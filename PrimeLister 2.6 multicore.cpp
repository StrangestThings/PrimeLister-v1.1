#include <algorithm>   // min, max
#include <cmath>       // sqrt, floor, log
#include <cstdint>     // uint8_t, uint64_t
#include <cstdlib>     // _dupenv_s, _putenv_s, free, atoi
#include <fstream>     // ofstream
#include <iomanip>     // setprecision, fixed
#include <iostream>    // cin, cout, cerr
#include <string>      // string
#include <vector>      // std::vector
#include <charconv>    // to_chars
#include <cstring>     // memcpy

// ---- OpenMP optional (Fallback without /openmp) ----
#ifdef _OPENMP
#include <omp.h>
#else
#include <chrono>
inline double omp_get_wtime() {
    static const auto t0 = std::chrono::steady_clock::now();
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}
inline int  omp_get_num_procs() { return 1; }
inline int  omp_get_max_threads() { return 1; }
inline void omp_set_dynamic(int) {}
inline void omp_set_num_threads(int) {}
inline int  omp_get_thread_num() { return 0; }
#endif

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <intrin.h>
#endif

using u64 = unsigned long long;

// ---------- helpers ----------
static inline u64 ceil_div_u64(u64 a, u64 b) {
    return (a + b - 1) / b;
}

static inline unsigned ctz64(uint64_t x) {
#if defined(_MSC_VER)
    // MSVC: use _BitScanForward64 (works without BMI1)
    unsigned long idx;
    _BitScanForward64(&idx, x);
    return (unsigned)idx;
#else
    return (unsigned)__builtin_ctzll(x);
#endif
}

// ---------- safe Env getter (avoids C4996) ----------
static bool get_env_str(const char* name, std::string& out) {
#if defined(_WIN32)
    char* buf = nullptr;
    size_t len = 0;
    if (_dupenv_s(&buf, &len, name) != 0 || !buf) return false;
    out.assign(buf, (len && buf[len - 1] == '\0') ? len - 1 : len);
    free(buf);
    return true;
#else
    const char* v = std::getenv(name);
    if (!v) return false;
    out.assign(v);
    return true;
#endif
}

// ---------- Thread setup: keep 'leave_free' threads free ----------
static int configure_threads_leave_free(int leave_free = 2) {
#ifdef _OPENMP
    int procs = omp_get_num_procs();
    int want = procs - leave_free;
    if (want < 1) want = 1;

    omp_set_dynamic(0);
    omp_set_num_threads(want);

#if defined(_WIN32)
    _putenv_s("OMP_PROC_BIND", "TRUE");
    _putenv_s("OMP_PLACES", "CORES");
#endif
    return want;
#else
    (void)leave_free;
    return 1;
#endif
}

// ---------- Fast output (buffered + to_chars); fixes remainder for cols>1 ----------
static inline void buf_append(std::string& buf, const char* s) {
    buf.append(s);
}
static inline void buf_append_char(std::string& buf, char c) {
    buf.push_back(c);
}
static inline void buf_append_u64(std::string& buf, u64 x) {
    char tmp[32];
    auto [ptr, ec] = std::to_chars(tmp, tmp + 32, x);
    (void)ec;
    buf.append(tmp, ptr);
}
static inline void buf_flush(std::ofstream& out, std::string& buf) {
    if (!buf.empty()) {
        out.write(buf.data(), (std::streamsize)buf.size());
        buf.clear();
    }
}

static void write_prime_table(const std::vector<u64>& P, u64 cols) {
    std::ofstream out("PrimeList.txt", std::ios::out | std::ios::binary);
    if (!out) { std::cerr << "Error: cannot open PrimeList.txt\n"; return; }

    if (cols < 1) cols = 1;

    // 8 MiB output buffer
    std::string buf;
    buf.reserve(8u * 1024u * 1024u);

    if (cols == 1) {
        // format: "index:\tprime\n"
        for (size_t i = 0; i < P.size(); ++i) {
            buf_append_u64(buf, (u64)(i + 1));
            buf_append(buf, ":\t");
            buf_append_u64(buf, P[i]);
            buf_append_char(buf, '\n');

            if (buf.size() >= 8u * 1024u * 1024u) buf_flush(out, buf);
        }
        buf_flush(out, buf);
        return;
    }

    auto header = [&](u64 start_index_1based) {
        u64 end_index_1based = std::min<u64>(start_index_1based + cols - 1, (u64)P.size());
        buf_append_u64(buf, start_index_1based);
        buf_append_char(buf, '-');
        buf_append_u64(buf, end_index_1based);
        buf_append_char(buf, '\t');
        };

    u64 start = 1;
    u64 count = 0;
    header(start);

    for (size_t i = 0; i < P.size(); ++i) {
        buf_append_u64(buf, P[i]);
        buf_append_char(buf, '\t');

        ++count;
        if (count == cols) {
            buf_append_char(buf, '\n');
            count = 0;
            start = (u64)(i + 2);                 // next 1-based index
            if (start <= (u64)P.size()) header(start);
        }

        if (buf.size() >= 8u * 1024u * 1024u) buf_flush(out, buf);
    }

    // If last row not complete, still end line
    if (count != 0) buf_append_char(buf, '\n');
    buf_flush(out, buf);
}

// ---------- Optimized formulas-mode (segmented + bitset + presieve p=5 & p=7) ----------
static std::vector<u64> primes_formeln(u64 N, u64 cols) {
    // Index range ~ N/3 plus a small safety margin (integer math)
    // Using ceil((N+2)/3) is safe for p = 3*i + 4/5
    u64 bereich = (N + 2ull) / 3ull + (15ull * std::max<u64>(1, cols));
    if (bereich == 0) return {};

    int T = omp_get_max_threads();
    if (T < 1) T = 1;

    // target: ~6 segments per thread for better load distribution
    u64 target_segments = (u64)T * 6ull;
    if (target_segments == 0) target_segments = 1;
    u64 SEG_IDX = (bereich + target_segments - 1) / target_segments;

    // Bitset memory guideline: 2–16 MiB for ~good cache behavior
    const u64 MIN_SEG = 1ull << 24;  // ~2 MiB bitset
    const u64 MAX_SEG = 1ull << 27;  // ~16 MiB bitset
    if (SEG_IDX < MIN_SEG) SEG_IDX = MIN_SEG;
    if (SEG_IDX > MAX_SEG) SEG_IDX = MAX_SEG;
    if (SEG_IDX > bereich) SEG_IDX = bereich;

    u64 nsegs = (bereich + SEG_IDX - 1) / SEG_IDX;

    // Global bounds (still useful as clamps); computed from full range
    const long long M = (long long)bereich;

    long long end1 = (M - 11) / 7; if (end1 < 0) end1 = -1;

    long long end2;
    if (M <= 4) {
        end2 = -1;
    }
    else {
        long double disc = 16.0L + 12.0L * (long double)M;
        long double root = (-8.0L + std::sqrt(disc)) / 6.0L;
        end2 = (long long)std::floor(root);
    }

    long long end3;
    if (M <= 7) {
        end3 = -1;
    }
    else {
        long double disc = 16.0L + 12.0L * (long double)M;
        long double root = (-10.0L + std::sqrt(disc)) / 6.0L;
        end3 = (long long)std::floor(root);
    }

    // Segment results
    std::vector<std::vector<u64>> buckets((size_t)nsegs);

    // telemetry (avoid false sharing)
    struct alignas(64) Counter { unsigned v = 0; };
    std::vector<Counter> segs_per_thread((size_t)T);

#pragma omp parallel
    {
        const int tid = omp_get_thread_num();

        // Thread-local reusable bitset words for max segment size
        std::vector<uint64_t> mark_words((size_t)((SEG_IDX + 63ull) >> 6), 0ull);

        auto setbit = [&](size_t idx) {
            mark_words[idx >> 6] |= (1ull << (idx & 63));
            };
        auto getbit = [&](size_t idx) -> bool {
            return (mark_words[idx >> 6] >> (idx & 63)) & 1ull;
            };

#pragma omp for schedule(guided, 1)
        for (long long s = 0; s < (long long)nsegs; ++s) {
            u64 base = (u64)s * SEG_IDX;
            u64 end = std::min<u64>(bereich, base + SEG_IDX);
            size_t seg_len = (size_t)(end - base);
            size_t nwords = (seg_len + 63u) >> 6;

            // clear only used words
            std::fill(mark_words.begin(), mark_words.begin() + nwords, 0ull);

            // ---- Per-segment bounds (big speed win) ----
            long long Mend = (long long)end;

            long long end1_seg = (Mend - 11) / 7;
            if (end1_seg < 0) end1_seg = -1;
            if (end1_seg > end1) end1_seg = end1;

            long long end2_seg;
            if (Mend <= 4) {
                end2_seg = -1;
            }
            else {
                long double disc = 16.0L + 12.0L * (long double)Mend;
                long double root = (-8.0L + std::sqrt(disc)) / 6.0L;
                end2_seg = (long long)std::floor(root);
                if (end2_seg > end2) end2_seg = end2;
            }

            long long end3_seg;
            if (Mend <= 7) {
                end3_seg = -1;
            }
            else {
                long double disc = 16.0L + 12.0L * (long double)Mend;
                long double root = (-10.0L + std::sqrt(disc)) / 6.0L;
                end3_seg = (long long)std::floor(root);
                if (end3_seg > end3) end3_seg = end3;
            }

            // --- Presieve p=5: i ≡ 0 (mod 10) [even], i ≡ 7 (mod 10) [odd]
            {
                u64 r10 = base % 10ull;

                u64 first0 = (r10 == 0ull) ? base : (base + (10ull - r10));
                for (u64 i = first0; i < end; i += 10ull) {
                    if (base == 0 && i == 0) continue; // don't kill p=5 (i=0 -> 5)
                    setbit((size_t)(i - base));
                }

                u64 add7 = (r10 <= 7 ? 7 - r10 : 17 - r10);
                u64 first7 = base + add7;
                for (u64 i = first7; i < end; i += 10ull) {
                    setbit((size_t)(i - base));
                }
            }

            // --- Presieve p=7: even i ≡ 10 (mod 14), odd i ≡ 1 (mod 14)
            {
                u64 r14 = base % 14ull;

                u64 add10 = (r14 <= 10 ? 10 - r14 : 24 - r14);
                u64 first10 = base + add10;
                for (u64 i = first10; i < end; i += 14ull) {
                    setbit((size_t)(i - base));
                }

                u64 add1 = (r14 <= 1 ? 1 - r14 : 15 - r14);
                u64 first1 = base + add1;
                for (u64 i = first1; i < end; i += 14ull) {
                    if (base == 0 && i == 1) continue; // don't kill p=7 (i=1 -> 7)
                    setbit((size_t)(i - base));
                }
            }

            // ----- Phase 1: j = 7*i + 10 ; step = 6*i + 10 ; i = 0,2,4,...
            for (long long i = 0; i <= end1_seg; i += 2) {
                u64 iu = (u64)i;
                u64 step = 6ull * iu + 10ull;
                u64 j0 = 7ull * iu + 10ull;

                if (j0 < base) {
                    u64 d = base - j0;
                    j0 += ceil_div_u64(d, step) * step;
                }
                for (u64 j = j0; j < end; j += step) setbit((size_t)(j - base));
            }

            // ----- Phase 2: j = 3*i^2 + 8*i + 4 ; step = 6*i + 8 ; i = 1,3,5,...
            for (long long i = 1; i <= end2_seg; i += 2) {
                u64 iu = (u64)i;
                u64 step = 6ull * iu + 8ull;
                u64 j0 = 3ull * iu * iu + 8ull * iu + 4ull;

                if (j0 < base) {
                    u64 d = base - j0;
                    j0 += ceil_div_u64(d, step) * step;
                }
                for (u64 j = j0; j < end; j += step) setbit((size_t)(j - base));
            }

            // ----- Phase 3: idx = 3*i^2 + 10*i + 7 ; step = 6*i + 10 ; i = 0,2,4,...
            for (long long i = 0; i <= end3_seg; i += 2) {
                u64 iu = (u64)i;
                u64 idx = 3ull * iu * iu + 10ull * iu + 7ull;
                u64 step = 6ull * iu + 10ull;

                if (idx < base) {
                    u64 d = base - idx;
                    u64 j0 = idx + ceil_div_u64(d, step) * step;
                    for (u64 j = j0; j < end; j += step) setbit((size_t)(j - base));
                }
                else if (idx < end) {
                    size_t loc = (size_t)(idx - base);
                    if (!getbit(loc)) setbit(loc);
                    for (u64 j = idx + step; j < end; j += step) setbit((size_t)(j - base));
                }
            }

            // ---- Collect candidates (word-wise with ctz) ----
            std::vector<u64> local;
            if (seg_len > 0) {
                // rough reserve (cheap heuristic)
                double approx = (double)seg_len / std::max(1.0, std::log((double)std::max<u64>(3, end)));
                local.reserve((size_t)approx + 32);
            }

            for (size_t wi = 0; wi < nwords; ++wi) {
                uint64_t inv = ~mark_words[wi];     // 1 bits = "unmarked"
                while (inv) {
                    unsigned b = ctz64(inv);
                    size_t off = (wi << 6) + (size_t)b;
                    if (off >= seg_len) break;

                    u64 iidx = base + (u64)off;
                    u64 p = ((iidx & 1ull) == 0ull) ? (3ull * iidx + 5ull) : (3ull * iidx + 4ull);
                    if (p <= N) local.push_back(p);

                    inv &= (inv - 1ull); // clear lowest set bit
                }
            }

            buckets[(size_t)s].swap(local);

            if ((size_t)tid < segs_per_thread.size()) segs_per_thread[(size_t)tid].v++;
        }
    } // omp parallel

    // ---- Merge (segment order), prepend 2 and 3 ----
    std::vector<u64> primes;
    size_t total = 1 + (N >= 3 ? 1 : 0);
    for (auto& b : buckets) total += b.size();
    primes.reserve(total);

    primes.push_back(2);
    if (N >= 3) primes.push_back(3);
    for (auto& b : buckets) primes.insert(primes.end(), b.begin(), b.end());

    // Telemetry
    std::cout << "\n[Telemetry] logical cores: " << omp_get_num_procs()
        << ", threads used: " << omp_get_max_threads()
        << ", segments: " << nsegs
        << ", seg_len≈ " << SEG_IDX << " indices (" << (SEG_IDX / 8) << " bytes bitset)\n";
    std::cout << "[Telemetry] segments per thread: ";
    for (size_t k = 0; k < segs_per_thread.size(); ++k) {
        std::cout << segs_per_thread[k].v << (k + 1 < segs_per_thread.size() ? ',' : '\n');
    }

    return primes;
}

// ------------------- MAIN -------------------
int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    std::cout << "****************************\n";
    std::cout << "  PrimeLister v2.6 (2026)   *\n";
    std::cout << "******************************\n\n";

    u64 N; u64 cols;
    std::cout << "Upper limit N: ";
    if (!(std::cin >> N)) return 0;
    std::cout << "Number of columns for PrimeList.txt: ";
    if (!(std::cin >> cols) || cols < 1) cols = 1;

    // --- Thread config: keep 2 free (override via ENV PRIMELISTER_THREADS) ---
    int T = 0;
#ifdef _OPENMP
    std::string env_threads;
    if (get_env_str("PRIMELISTER_THREADS", env_threads)) {
        int v = std::atoi(env_threads.c_str());
        if (v > 0) { omp_set_dynamic(0); omp_set_num_threads(v); T = v; }
    }
#endif
    if (T == 0) T = configure_threads_leave_free(2);

#ifdef _OPENMP
    std::cout << "Detected logical cores: " << omp_get_num_procs()
        << " | using threads: " << T << "\n";
#else
    std::cout << "OpenMP not enabled -> running single-threaded.\n";
#endif

    double t0 = omp_get_wtime();
    std::vector<u64> P = primes_formeln(N, cols);
    double t1 = omp_get_wtime();

    // Short console output
    std::cout << "\nFirst 100 primes:\n";
    for (size_t i = 0; i < std::min<size_t>(100, P.size()); ++i) std::cout << P[i] << ' ';

    std::cout << "\n\nLast 100 primes:\n";
    for (size_t i = (P.size() > 100 ? P.size() - 100 : 0); i < P.size(); ++i) std::cout << P[i] << ' ';

    std::cout << "\n\nOverall primes found: " << P.size() << "\n";
    std::cout << "Overall calculation time [s]: " << std::fixed << std::setprecision(3) << (t1 - t0) << "\n";
    if (t1 - t0 > 0) {
        auto rate = (double)P.size() / (t1 - t0);
        std::cout << "~ " << (u64)rate << " primes/s\n";
    }

    std::cout << "\nWriting PrimeList.txt ...\n";
    write_prime_table(P, cols);
    std::cout << "Done. Press Enter to exit...";
    std::cin.get(); std::cin.get();
    return 0;
}

