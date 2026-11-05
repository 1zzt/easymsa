#define MODE 3
#include "internal.hpp"
#include "msakit/core.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string_view>
#include <unordered_map>
#include <vector>

// MODE is fixed per project: 1 fast-dedup, 2 kmer-center, 3 affine-star.
namespace msakit {
namespace {
using Score = std::int64_t;
struct Projection { std::string at_center; std::vector<std::string> insertions; };
bool residue(char c) { return std::string_view("ACGTUN").find(c) != std::string_view::npos; }
Score sub(char a, char b) {
    if (a == 'N' || b == 'N') return 0;
    if (a == 'U') a = 'T';
    if (b == 'U') b = 'T';
    return a == b ? (MODE == 3 ? 2 : 1) : -1;
}
Projection project(std::string a, std::string b, std::size_t length) {
    std::reverse(a.begin(), a.end());
    std::reverse(b.begin(), b.end());
    Projection p{std::string(length, '-'), std::vector<std::string>(detail::checked_add(length, 1))};
    std::size_t pos = 0;
    for (std::size_t k = 0; k < a.size(); ++k) {
        if (a[k] == '-') p.insertions[pos].push_back(b[k]);
        else p.at_center[pos++] = b[k];
    }
    return p;
}
Projection pair_linear(std::string_view a, std::string_view b) {
    detail::check_score_bounds(a.size(), b.size());
    const auto h = detail::checked_add(a.size(), 1), w = detail::checked_add(b.size(), 1);
    const auto cells = detail::checked_multiply(h, w);
    // One byte per traceback cell and two score rows; preserve upstream ties.
    std::vector<std::uint8_t> trace(cells, 0);
    std::vector<Score> prev(w), curr(w);
    for (std::size_t j = 1; j < w; ++j) { prev[j] = -static_cast<Score>(j); trace[j] = 2; }
    for (std::size_t i = 1; i < h; ++i) {
        curr[0] = -static_cast<Score>(i); trace[i*w] = 1;
        for (std::size_t j = 1; j < w; ++j) {
            auto best = prev[j-1] + sub(a[i-1], b[j-1]);
            std::uint8_t direction = 0;
            if (prev[j]-1 > best) { best = prev[j]-1; direction = 1; }
            if (curr[j-1]-1 > best) { best = curr[j-1]-1; direction = 2; }
            curr[j] = best; trace[i*w+j] = direction;
        }
        prev.swap(curr);
    }
    std::string x, y;
    x.reserve(detail::checked_add(a.size(), b.size())); y.reserve(x.capacity());
    auto i = a.size(), j = b.size();
    while (i || j) {
        const auto d = trace[i*w+j];
        if (d == 0) { x.push_back(a[--i]); y.push_back(b[--j]); }
        else if (d == 1) { x.push_back(a[--i]); y.push_back('-'); }
        else { x.push_back('-'); y.push_back(b[--j]); }
    }
    return project(std::move(x), std::move(y), a.size());
}
Projection pair_affine(std::string_view a, std::string_view b) {
    const auto total = detail::checked_add(a.size(), b.size());
    if (total > static_cast<std::size_t>(std::numeric_limits<Score>::max()/8))
        throw detail::Error(detail::kCoreError, "affine score range exceeded");
    const auto h = detail::checked_add(a.size(), 1), w = detail::checked_add(b.size(), 1);
    const auto cells = detail::checked_multiply(h, w);
    (void)detail::checked_multiply(cells, 3);
    constexpr Score neg = std::numeric_limits<Score>::min()/4;
    // State 0 consumes both; 1 consumes center; 2 consumes query.
    // Gap length k costs 3 + (k-1), including terminal gaps.
    std::array<std::vector<Score>,3> p, c;
    std::array<std::vector<std::uint8_t>,3> trace;
    for (int s=0;s<3;++s) { p[s].assign(w,neg); c[s].assign(w,neg); trace[s].assign(cells,0); }
    p[0][0]=0;
    for (std::size_t j=1;j<w;++j) { p[2][j]=-3-static_cast<Score>(j-1); trace[2][j]=j==1?0:2; }
    const auto choose = [](const std::array<Score,3>& v) {
        return static_cast<std::uint8_t>(std::max_element(v.begin(), v.end())-v.begin());
    };
    for (std::size_t i=1;i<h;++i) {
        c[0][0]=neg; c[1][0]=-3-static_cast<Score>(i-1); c[2][0]=neg;
        trace[1][i*w]=i==1?0:1;
        for (std::size_t j=1;j<w;++j) {
            std::array<Score,3> v{p[0][j-1],p[1][j-1],p[2][j-1]};
            auto s=choose(v); c[0][j]=v[s]+sub(a[i-1],b[j-1]); trace[0][i*w+j]=s;
            v={p[0][j]-3,p[1][j]-1,p[2][j]-3};
            s=choose(v); c[1][j]=v[s]; trace[1][i*w+j]=s;
            v={c[0][j-1]-3,c[1][j-1]-3,c[2][j-1]-1};
            s=choose(v); c[2][j]=v[s]; trace[2][i*w+j]=s;
        }
        p.swap(c);
    }
    std::uint8_t s=choose({p[0].back(),p[1].back(),p[2].back()});
    std::string x,y; x.reserve(total); y.reserve(total);
    auto i=a.size(),j=b.size();
    while (i||j) {
        const auto next=trace[s][i*w+j];
        if (s==0) { x.push_back(a[--i]); y.push_back(b[--j]); }
        else if (s==1) { x.push_back(a[--i]); y.push_back('-'); }
        else { x.push_back('-'); y.push_back(b[--j]); }
        s=next;
    }
    return project(std::move(x),std::move(y),a.size());
}
std::size_t choose_center(const std::vector<SequenceRecord>& input) {
    if (MODE != 2) return 0;
    // Linear-time heuristic: normalized 4-mer profile nearest the mean.
    // This is not an exact pairwise-distance medoid.
    std::vector<std::array<double,256>> profiles(input.size());
    std::array<double,256> mean{};
    double mean_length=0;
    for (std::size_t i=0;i<input.size();++i) {
        auto& f=profiles[i]; f.fill(0);
        unsigned code=0,valid=0; std::size_t count=0;
        for (char ch:input[i].sequence) {
            unsigned b=0;
            if (ch=='C') b=1; else if(ch=='G') b=2; else if(ch=='T'||ch=='U') b=3;
            else if(ch!='A') {code=0;valid=0;continue;}
            code=((code<<2)|b)&255; ++valid;
            if(valid>=4) {f[code]+=1;++count;}
        }
        for(std::size_t k=0;k<256;++k) { if(count) f[k]/=static_cast<double>(count); mean[k]+=f[k]/input.size(); }
        mean_length+=static_cast<double>(input[i].sequence.size())/input.size();
    }
    std::size_t best=0; double best_distance=std::numeric_limits<double>::infinity();
    for(std::size_t i=0;i<input.size();++i) {
        const auto delta=(static_cast<double>(input[i].sequence.size())-mean_length)/std::max(1.0,mean_length);
        double distance=0.05*delta*delta;
        for(std::size_t k=0;k<256;++k) { auto d=profiles[i][k]-mean[k]; distance+=d*d; }
        if(distance<best_distance || (distance==best_distance && input[i].sequence<input[best].sequence)) {
            best=i;best_distance=distance;
        }
    }
    return best;
}
} // namespace
Alignment align(const std::vector<SequenceRecord>& input) {
    if(input.empty()) throw detail::Error(detail::kCoreError,"at least one record is required");
    for(const auto& r:input) if(r.sequence.empty() || !std::all_of(r.sequence.begin(),r.sequence.end(),residue))
        throw detail::Error(detail::kCoreError,"nonempty cleaned ACGTUN sequences are required");
    if(input.size()==1) return Alignment{{input[0].sequence}};
    const auto center_index=choose_center(input);
    const auto& center=input[center_index].sequence;
    const auto slots=detail::checked_add(center.size(),1);
    std::vector<std::size_t> widths(slots,0), mapping(input.size());
    std::unordered_map<std::string_view,std::size_t> seen;
    std::vector<Projection> projections;
    // Cache by complete cleaned sequence, never by FASTA header.
    for(std::size_t row=0;row<input.size();++row) {
        const auto& seq=input[row].sequence;
        auto found=seen.find(seq);
        if(found!=seen.end()) {mapping[row]=found->second;continue;}
        Projection p;
        if(seq==center) p={center,std::vector<std::string>(slots)};
        else if(MODE==3) p=pair_affine(center,seq);
        else p=pair_linear(center,seq);
        mapping[row]=projections.size(); seen.emplace(seq,projections.size());
        for(std::size_t s=0;s<slots;++s) widths[s]=std::max(widths[s],p.insertions[s].size());
        projections.push_back(std::move(p));
    }
    auto width=center.size();
    for(auto n:widths) width=detail::checked_add(width,n);
    (void)detail::checked_multiply(input.size(),width);
    Alignment result; result.rows.resize(input.size());
    for(std::size_t row=0;row<input.size();++row) {
        const auto& p=projections[mapping[row]]; auto& out=result.rows[row]; out.reserve(width);
        for(std::size_t s=0;s<slots;++s) {
            out+=p.insertions[s]; out.append(widths[s]-p.insertions[s].size(),'-');
            if(s<center.size()) out.push_back(p.at_center[s]);
        }
    }
    return result;
}
} // namespace msakit
