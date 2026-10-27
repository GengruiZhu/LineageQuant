// LineageQuant 1.1.0 -- reference preparation (prep.cpp), shared with the pipeline
#pragma once
#include "binio.hpp"
#include "common.hpp"

#include <functional>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

namespace lq {

std::string_view py_strip(std::string_view s);
std::string_view py_rstrip(std::string_view s);
std::vector<std::string_view> py_split(std::string_view s, char sep);
std::vector<std::string_view> py_split_ws(std::string_view s);
std::string_view before_dot(std::string_view s);
bool py_int(std::string_view s, long long& v);
void remove_dir_tree(const std::string& d);
std::string make_temp_dir(const std::string& base, const std::string& prefix);
void concat_files(const std::vector<std::string>& in, const std::string& out);
bool is_low_complexity(const std::string& km, int k, int thr);

// build_markers.build_markers(): returns {lineage: OUTDIR/lineage.fasta}; *order = lineages in file-creation order
std::map<std::string, std::string> build_markers(const std::vector<std::pair<std::string, std::string>>& anc,
                                                 const std::string& outdir, int k, int threads, double at_gc_frac,
                                                 const std::vector<std::pair<std::string, std::string>>& poly,
                                                 const std::string& mode, bool keep_intermediate,
                                                 std::vector<std::string>* order);
// GenomeAnnotationPair(fasta, gff, workdir, name).prepare().extract_transcripts(out_fa, id_prefix)
void genome_pair_extract(const std::string& fasta, const std::string& gff, const std::string& workdir, const std::string& name,
                         const std::string& out_fa, const std::string& id_prefix, long long min_chr_len = 1000000);
bool write_fai(const std::string& fa, const std::string& fai, std::vector<std::pair<std::string, long long>>* lens, std::string* err);

}  // namespace lq
