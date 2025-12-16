#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <filesystem>
#include <algorithm>
#include <cstdint>

using namespace std;
namespace fs = std::filesystem;
using u64 = uint64_t;

class HyperLogLog {
public:
    explicit HyperLogLog(uint8_t p_) : p(p_), m(1u << p_), M(m, 0) {
        if (m == 16)      alpha = 0.673;
        else if (m == 32) alpha = 0.697;
        else if (m == 64) alpha = 0.709;
        else              alpha = 0.7213 / (1.0 + 1.079 / m);
    }

    void add(u64 h) {
        if (h == 0) h = 1;
        uint32_t idx = h >> (64 - p);
        u64 w = h << p;
        int lz = __builtin_clzll(w);
        int rho = lz + 1;
        if (rho > M[idx]) {
            M[idx] = (uint8_t)rho;
        }
    }

    double estimate() const {
        double sum = 0.0;
        int V = 0;
        for (uint32_t j = 0; j < m; ++j) {
            sum += std::ldexp(1.0, -M[j]);
            if (M[j] == 0) ++V;
        }
        double E = alpha * (double)m * (double)m / sum;
        if (E <= 5.0 * m && V > 0) {
            E = m * std::log((double)m / V);
        }
        return E;
    }

    bool save(const fs::path& out_path) const {
        ofstream out(out_path, ios::binary);
        if (!out) {
            cerr << "[WARN] No pude crear sketch: " << out_path << "\n";
            return false;
        }
        const char magic[4] = {'H','L','L','1'};
        out.write(magic, 4);
        out.write(reinterpret_cast<const char*>(&p), sizeof(p));
        out.write(reinterpret_cast<const char*>(&m), sizeof(m));
        out.write(reinterpret_cast<const char*>(M.data()), M.size());
        if (!out.good()) {
            cerr << "[WARN] Error al escribir sketch: " << out_path << "\n";
            return false;
        }
        return true;
    }

    static HyperLogLog load(const fs::path& in_path) {
        ifstream in(in_path, ios::binary);
        if (!in) {
            throw runtime_error("No puedo abrir sketch: " + in_path.string());
        }
        char magic[4];
        in.read(magic, 4);
        if (in.gcount() != 4 || magic[0] != 'H' || magic[1] != 'L' ||
            magic[2] != 'L' || magic[3] != '1') {
            throw runtime_error("Formato de sketch inválido en: " + in_path.string());
        }
        uint8_t p_read;
        uint32_t m_read;
        in.read(reinterpret_cast<char*>(&p_read), sizeof(p_read));
        in.read(reinterpret_cast<char*>(&m_read), sizeof(m_read));
        HyperLogLog hll(p_read);
        if (hll.m != m_read) {
            throw runtime_error("m inconsistente al cargar sketch: " + in_path.string());
        }
        in.read(reinterpret_cast<char*>(hll.M.data()), hll.M.size());
        if (!in.good()) {
            throw runtime_error("Error leyendo registros M en: " + in_path.string());
        }
        return hll;
    }

    uint8_t getP() const { return p; }
    uint32_t getM() const { return m; }
    const vector<uint8_t>& getRegisters() const { return M; }

private:
    uint8_t p;
    uint32_t m;
    double alpha;
    vector<uint8_t> M;
};

void process_minimizer_file(const fs::path& file_path, HyperLogLog& hll) {
    ifstream in(file_path);
    if (!in) {
        cerr << "[WARN] No pude abrir minimizers: " << file_path << "\n";
        return;
    }
    cout << "[INFO] Leyendo " << file_path.filename() << "\n";
    u64 h;
    size_t idx;
    while (in >> h >> idx) {
        hll.add(h);
    }
}

int main() {
    string input_dir_str, output_dir_str;
    string pattern_suffix = ".txt";

    cout << "Carpeta con archivos de minimizers (*.txt): ";
    cin >> input_dir_str;

    cout << "Carpeta de salida para sketches (.hll): ";
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

    const uint8_t P = 14;

    for (const auto& entry : fs::directory_iterator(input_dir)) {
        if (!entry.is_regular_file()) continue;

        fs::path pth = entry.path();
        string fname = pth.filename().string();

        if (fname.size() < pattern_suffix.size()) continue;
        if (fname.substr(fname.size() - pattern_suffix.size()) != pattern_suffix)
            continue;

        HyperLogLog hll(P);
        process_minimizer_file(pth, hll);

        double N_est = hll.estimate();
        cout << "[RESULT] " << fname
             << "  cardinalidad_estimada=" << N_est << "\n";

        string stem = pth.stem().string();
        fs::path out_sketch = output_dir / (stem + ".hll");

        if (hll.save(out_sketch)) {
            cout << "    -> Sketch guardado en: " << out_sketch.filename() << "\n";
            try {
                uintmax_t size_bytes = fs::file_size(out_sketch);
                double size_kib = static_cast<double>(size_bytes) / 1024.0;
                uint32_t m_reg = hll.getM();
                uint64_t theoretical_bytes = 4 + 1 + 4 + m_reg;
                cout << "       Tamaño archivo: "
                     << size_bytes << " bytes ("
                     << size_kib << " KiB)\n";
                cout << "       Tamaño teorico: "
                     << theoretical_bytes << " bytes"
                     << "  [m = " << m_reg << "]\n";
            } catch (const std::exception& e) {
                cerr << "       [WARN] No pude obtener file_size: " << e.what() << "\n";
            }
        }
    }

    cout << "[OK] Construcción de sketches completada.\n";
    return 0;
}
