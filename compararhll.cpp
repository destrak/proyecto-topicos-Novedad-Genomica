#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <filesystem>
#include <stdexcept>
#include <unordered_set>
#include <string>

using namespace std;
namespace fs = std::filesystem;

class HyperLogLog {
public:
    using u64 = uint64_t;

    explicit HyperLogLog(uint8_t p_) : p(p_), m(1u << p_), M(m, 0) {}

    static HyperLogLog load(const fs::path& in_path) {
        ifstream in(in_path, ios::binary);
        if (!in)
            throw runtime_error("No puedo abrir: " + in_path.string());

        char magic[4];
        in.read(magic, 4);
        if (in.gcount() != 4 ||
            magic[0] != 'H' || magic[1] != 'L' ||
            magic[2] != 'L' || magic[3] != '1')
            throw runtime_error("Formato inválido en: " + in_path.string());

        uint8_t p_read;
        uint32_t m_read;
        in.read((char*)&p_read, sizeof(p_read));
        in.read((char*)&m_read, sizeof(m_read));

        HyperLogLog hll(p_read);
        if (hll.m != m_read)
            throw runtime_error("m inconsistente en: " + in_path.string());

        in.read((char*)hll.M.data(), hll.M.size());
        if (!in.good())
            throw runtime_error("Error leyendo M[] en: " + in_path.string());

        return hll;
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

        if (E <= 5.0 * m && V > 0)
            E = m * log((double)m / V);

        return E;
    }

    static HyperLogLog make_union(const HyperLogLog& R, const HyperLogLog& S) {
        if (R.p != S.p || R.m != S.m)
            throw runtime_error("Sketches incompatibles (p/m).");

        HyperLogLog U(R.p);
        for (uint32_t j = 0; j < R.m; ++j)
            U.M[j] = std::max(R.M[j], S.M[j]);
        return U;
    }

    uint8_t p;
    uint32_t m;
    vector<uint8_t> M;
};

int main(int argc, char** argv) {
    if (argc != 5) {
        cerr << "Uso: " << argv[0]
             << " sketch_R.hll carpeta_S usados.txt result.csv\n";
        return 1;
    }

    fs::path pathR   = argv[1];
    fs::path dirS    = argv[2];
    fs::path usedLst = argv[3];
    fs::path outCSV  = argv[4];

    if (!fs::exists(pathR) || !fs::is_regular_file(pathR)) {
        cerr << "Error: sketch_R no existe o no es archivo regular.\n";
        return 1;
    }
    if (!fs::exists(dirS) || !fs::is_directory(dirS)) {
        cerr << "Error: carpeta_S no existe o no es directorio.\n";
        return 1;
    }
    if (!fs::exists(usedLst) || !fs::is_regular_file(usedLst)) {
        cerr << "Error: usados.txt no existe o no es archivo regular.\n";
        return 1;
    }

    try {
        HyperLogLog R = HyperLogLog::load(pathR);
        double sizeR  = R.estimate();

        cout << "[INFO] Catalogo R: " << pathR << "\n";
        cout << "[INFO] |R| ≈ " << sizeR << "\n";

        unordered_set<string> used_names;
        {
            ifstream in_used(usedLst);
            if (!in_used) {
                cerr << "Error: no pude abrir lista de usados: " << usedLst << "\n";
                return 1;
            }
            string line;
            while (getline(in_used, line)) {
                if (!line.empty())
                    used_names.insert(line);
            }
        }

        ofstream out(outCSV);
        if (!out) {
            cerr << "Error: no pude crear " << outCSV << "\n";
            return 1;
        }

        out << "sketch_S,|R|,|S|,|R_union_S|,|S_minus_R|,rho\n";

        int count_outside = 0;

        for (const auto& entry : fs::directory_iterator(dirS)) {
            if (!entry.is_regular_file()) continue;

            fs::path pS = entry.path();
            if (pS.extension() != ".hll") continue;

            string fname = pS.filename().string();

            if (fs::equivalent(pS, pathR))
                continue;

            if (used_names.find(fname) != used_names.end())
                continue;

            try {
                HyperLogLog S = HyperLogLog::load(pS);
                double sizeS  = S.estimate();

                HyperLogLog U = HyperLogLog::make_union(R, S);
                double sizeU  = U.estimate();

                double sizeS_minus_R = sizeU - sizeR;
                if (sizeS_minus_R < 0) sizeS_minus_R = 0.0;

                double rho = (sizeS > 0.0) ? (sizeS_minus_R / sizeS) : 0.0;

                cout << "\n[GENOMA FUERA DEL CATALOGO] " << fname << "\n";
                cout << "  |S|      ≈ " << sizeS << "\n";
                cout << "  |R ∪ S|  ≈ " << sizeU << "\n";
                cout << "  |S\\R|    ≈ " << sizeS_minus_R << "\n";
                cout << "  rho      ≈ " << rho << "\n";

                out << fname << ","
                    << sizeR << ","
                    << sizeS << ","
                    << sizeU << ","
                    << sizeS_minus_R << ","
                    << rho << "\n";

                ++count_outside;

            } catch (const std::exception& exS) {
                cerr << "[WARN] Error procesando " << pS
                     << ": " << exS.what() << "\n";
            }
        }

        cout << "\n[OK] Resultados escritos en: " << outCSV << "\n";
        cout << "[INFO] Genomas fuera del catalogo procesados: "
             << count_outside << "\n";

    } catch (const std::exception& ex) {
        cerr << "Error global: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}
