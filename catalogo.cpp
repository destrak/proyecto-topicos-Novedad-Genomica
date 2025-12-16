#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <filesystem>
#include <algorithm>
#include <random>
#include <stdexcept>
#include <cstdint>
#include <cmath>

using namespace std;
namespace fs = std::filesystem;

class HyperLogLog {
public:
    using u64 = uint64_t;

    explicit HyperLogLog(uint8_t p_) : p(p_), m(1u << p_), M(m, 0) {}

    static HyperLogLog load(const fs::path& in_path) {
        ifstream in(in_path, ios::binary);
        if (!in) {
            throw runtime_error("No puedo abrir sketch: " + in_path.string());
        }

        char magic[4];
        in.read(magic, 4);
        if (in.gcount() != 4 ||
            magic[0] != 'H' || magic[1] != 'L' ||
            magic[2] != 'L' || magic[3] != '1') {
            throw runtime_error("Formato inválido (magic) en: " + in_path.string());
        }

        uint8_t p_read;
        uint32_t m_read;
        in.read(reinterpret_cast<char*>(&p_read), sizeof(p_read));
        in.read(reinterpret_cast<char*>(&m_read), sizeof(m_read));

        HyperLogLog hll(p_read);
        if (hll.m != m_read) {
            throw runtime_error("m inconsistente al cargar: " + in_path.string());
        }

        in.read(reinterpret_cast<char*>(hll.M.data()), hll.M.size());
        if (!in.good()) {
            throw runtime_error("Error leyendo M[] en: " + in_path.string());
        }

        return hll;
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

    double estimate() const {
        double sum = 0.0;
        int V = 0;

        for (uint8_t x : M) {
            sum += std::ldexp(1.0, -x);
            if (x == 0) ++V;
        }

        double alpha;
        if (m == 16)      alpha = 0.673;
        else if (m == 32) alpha = 0.697;
        else if (m == 64) alpha = 0.709;
        else              alpha = 0.7213 / (1.0 + 1.079 / m);

        double E = alpha * (double)m * (double)m / sum;

        if (E <= 5.0 * m && V > 0) {
            E = m * log((double)m / V);
        }

        return E;
    }

    void merge(const HyperLogLog& other) {
        if (p != other.p || m != other.m) {
            throw runtime_error("No se puede hacer merge: p/m distintos.");
        }
        for (uint32_t i = 0; i < m; ++i) {
            if (other.M[i] > M[i]) {
                M[i] = other.M[i];
            }
        }
    }

    uint8_t p;
    uint32_t m;
    vector<uint8_t> M;
};

int main() {
    string input_dir_str;
    string output_file_str;
    int N;

    cout << "Carpeta con archivos HLL (*.hll): ";
    cin >> input_dir_str;

    cout << "Archivo de salida para el catalogo (ej. catalogo_7001_genomas.hll): ";
    cin >> output_file_str;

    cout << "Numero de archivos a seleccionar para el catalogo (N): ";
    cin >> N;

    fs::path input_dir(input_dir_str);
    if (!fs::exists(input_dir) || !fs::is_directory(input_dir)) {
        cerr << "Error: carpeta de entrada invalida.\n";
        return 1;
    }

    const string suffix = ".hll";
    vector<fs::path> hll_files;

    for (const auto& entry : fs::directory_iterator(input_dir)) {
        if (!entry.is_regular_file()) continue;
        fs::path p = entry.path();
        string fname = p.filename().string();
        if (fname.size() >= suffix.size() &&
            fname.substr(fname.size() - suffix.size()) == suffix) {
            hll_files.push_back(p);
        }
    }

    if (hll_files.empty()) {
        cerr << "No se encontraron archivos *.hll en la carpeta.\n";
        return 1;
    }

    cout << "[INFO] Se encontraron " << hll_files.size()
         << " archivos HLL.\n";

    random_device rd;
    mt19937 rng(rd());
    shuffle(hll_files.begin(), hll_files.end(), rng);

    if ((int)hll_files.size() < N) {
        cout << "[WARN] Solo hay " << hll_files.size()
             << " archivos. Se usaran todos para el catalogo.\n";
        N = (int)hll_files.size();
    } else {
        cout << "[INFO] Se seleccionaran " << N
             << " archivos aleatorios para el catalogo.\n";
    }

    vector<fs::path> chosen(hll_files.begin(), hll_files.begin() + N);
    vector<fs::path> not_chosen(hll_files.begin() + N, hll_files.end());

    string usados_list = output_file_str + "_usados.txt";
    string no_usados_list = output_file_str + "_no_usados.txt";

    ofstream out_usados(usados_list);
    ofstream out_no_usados(no_usados_list);

    if (!out_usados || !out_no_usados) {
        cerr << "Error: no pude crear los archivos de usados/no usados.\n";
        return 1;
    }

    cout << "[INFO] Cargando y uniendo sketches del catalogo...\n";

    HyperLogLog catalogo = HyperLogLog::load(chosen[0]);
    out_usados << chosen[0].filename().string() << "\n";

    for (size_t i = 1; i < chosen.size(); ++i) {
        out_usados << chosen[i].filename().string() << "\n";
        HyperLogLog h = HyperLogLog::load(chosen[i]);
        catalogo.merge(h);
    }

    for (const auto& p : not_chosen) {
        out_no_usados << p.filename().string() << "\n";
    }

    catalogo.save(output_file_str);

    cout << "\n[OK] Catalogo construido con " << N << " genomas\n";
    cout << "    Archivo: " << output_file_str << "\n";
    cout << "    Estimacion |R| ≈ " << catalogo.estimate() << "\n";
    cout << "    Lista usados:    " << usados_list << "\n";
    cout << "    Lista no usados: " << no_usados_list << "\n";

    return 0;
}
