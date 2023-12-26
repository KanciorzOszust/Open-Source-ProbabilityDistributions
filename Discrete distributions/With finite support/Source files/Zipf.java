public class Zipf extends MathLibrary{
    double s;
    int N;

    Zipf(double s, int N) {
        if (s < 0) throw new IllegalArgumentException("s >= 0");
        if (N < 1) throw new IllegalArgumentException("N >= 1, k >= 1");
        this.s = s;
        this.N = N;
    }

    double PMF(int k) {
        return 1.0 / doublePower(k, s) / generalizedHarmonicNumber(N, s);
    }

    double CDF(int k) {
        return generalizedHarmonicNumber(k, s) / generalizedHarmonicNumber(N, s);
    }

    double Mean() {
        return generalizedHarmonicNumber(N, s - 1) / generalizedHarmonicNumber(N, s);
    }

    int Mode() {
        return 1;
    }

    double Variance() {
        double part1 = generalizedHarmonicNumber(N, s - 2) / generalizedHarmonicNumber(N, s);
        double part2 = power(generalizedHarmonicNumber(N, s - 1), 2) / power(generalizedHarmonicNumber(N, s), 2);
        return part1 - part2;
    }

    double Entropy(int k) {
        double wartosc = s / generalizedHarmonicNumber(N, s);
        for (int i = 1; i < N; i++) {
            wartosc += ln(k) / doublePower(k, s);
        }
        return wartosc + ln(generalizedHarmonicNumber(N, s));
    }

    double MGF(double t) {
        double wartosc = 1.0 / generalizedHarmonicNumber(N, s);
        for (int i = 1; i < N; i++) {
            wartosc += euler(i * t) / doublePower(i, s);
        }
        return wartosc;
    }
}
