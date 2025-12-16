#include <bits/stdc++.h>
using namespace std;
namespace fs = std::filesystem;
using u64 = uint64_t;

inline int base2bits(char c) {
    switch (std::toupper(static_cast<unsigned char>(c))) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return -1;
    }
}

inline u64 hash64(u64 x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    x ^= (x >> 31);
    return x;
}

struct MinEntry {
    size_t kmer_idx;
    u64 h;
};

void process_sequence_minimizers(const string& seq,
                                 int k,
                                 int w,
                                 ofstream& out)
{
    if (k <= 0 || w <= 0) return;

    u64 mask;
    if (2 * k >= 64) mask = ~0ULL;
    else            mask = (1ULL << (2 * k)) - 1ULL;

    u64 fwd = 0;
    u64 rev = 0;
    int len = 0;
    size_t kmer_idx = 0;

    deque<MinEntry> dq;
    bool has_last_min = false;
    u64 last_min_hash = 0;

    for (size_t i = 0; i < seq.size(); ++i) {
        int code = base2bits(seq[i]);
        if (code < 0) {
            len = 0;
            fwd = 0;
            rev = 0;
            dq.clear();
            continue;
        }

        fwd = ((fwd << 2) | (u64)code) & mask;

        u64 comp = (u64)(3 - code);
        rev = (rev >> 2) | (comp << (2 * (k - 1)));

        if (len < k) {
            ++len;
            if (len < k) continue;
        }

        u64 canon = (fwd < rev) ? fwd : rev;
        u64 h = hash64(canon);

        size_t idx = kmer_idx++;

        size_t min_valid_idx = (idx >= (size_t)w - 1) ? idx - (w - 1) : 0;
        while (!dq.empty() && dq.front().kmer_idx < min_valid_idx) {
            dq.pop_front();
        }

        while (!dq.empty() && dq.back().h >= h) {
            dq.pop_back();
        }

        dq.push_back({idx, h});

        if (idx + 1 < (size_t)w) continue;

        u64 current_min = dq.front().h;

        if (!has_last_min || current_min != last_min_hash) {
            has_last_min = true;
            last_min_hash = current_min;
            out << current_min << '\t' << dq.front().kmer_idx << '\n';
        }
    }
}

void process_fasta_minimizers(const fs::path& fasta_path,
                              const fs::path& out_path,
                              int k,
                              int w)
{
    ifstream in(fasta_path);
    if (!in) {
        cerr << "[WARN] No pude abrir: " << fasta_path << "\n";
        return;
    }

    ofstream out(out_path);
    if (!out) {
        cerr << "[WARN] No pude crear: " << out_path << "\n";
        return;
    }

    cout << "[INFO] Procesando " << fasta_path.filename()
         << " -> " << out_path.filename() << "\n";

    string line;
    string current_seq;

    while (getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!current_seq.empty()) {
                process_sequence_minimizers(current_seq, k, w, out);
                current_seq.clear();
            }
        } else {
            for (char c : line) {
                if (!isspace((unsigned char)c))
                    current_seq.push_back(c);
            }
        }
    }
    if (!current_seq.empty()) {
        process_sequence_minimizers(current_seq, k, w, out);
    }
}

int main() {
    int k, w;
    string input_dir_str, output_dir_str;

    cout << "k (21 o 31): ";
    cin >> k;
    cout << "w (tamaño ventana para minimizers, ej. 10, 20, 50): ";
    cin >> w;

    cout << "Carpeta de entrada: ";
    cin >> input_dir_str;
    cout << "Carpeta de salida: ";
    cin >> output_dir_str;

    fs::path input_dir(input_dir_str);
    fs::path output_dir(output_dir_str);

    if (!fs::exists(input_dir) || !fs::is_directory(input_dir)) {
        cerr << "Error: carpeta de entrada inválida.\n";
        return 1;
    }
    if (!fs::exists(output_dir)) {
        fs::create_directories(output_dir);
    }

    for (const auto& entry : fs::directory_iterator(input_dir)) {
        if (!entry.is_regular_file()) continue;

        fs::path p = entry.path();
        string ext = p.extension().string();
        transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
        if (ext != ".fna" && ext != ".fa" && ext != ".fasta")
            continue;

        string stem = p.stem().string();
        fs::path out_file = output_dir /
            (stem + "_k" + to_string(k) + "_w" + to_string(w) + "_minimizers.txt");

        process_fasta_minimizers(p, out_file, k, w);
    }

    cout << "[OK] Minimizers generados.\n";
    return 0;
}
