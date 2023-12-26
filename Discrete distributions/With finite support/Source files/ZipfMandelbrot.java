public class ZipfMandelbrot extends MathLibrary{
    int N;
    double q;
    double s;

    ZipfMandelbrot(int N, double q, double s) {
        if (N < 1) throw new IllegalArgumentException("N >= 1");
        if (q < 0) throw new IllegalArgumentException("q >= 0");
        if (s <= 0) throw new IllegalArgumentException("s > 0");
        this.N = N;
        this.q = q;
        this.s = s;
    }

    double PMF(int k) {
        if (k < 1 || k > N) throw new IllegalArgumentException("1 <= k < = N");
        return 1 / doublePower(k + q, s) / parameterHarmonicNumber(N, q, s);
    }

    double CDF(int k) {
        if (k < 1 || k > N) throw new IllegalArgumentException("1 <= k < = N");
        return parameterHarmonicNumber(k, q, s) / parameterHarmonicNumber(N, q, s);
    }

    double Mean() {
        return parameterHarmonicNumber(N, q, s - 1) / parameterHarmonicNumber(N, q, s) - q;
    }

    int Mode() {
        return 1;
    }

    double Entropy() {
        double value = 0;
        for (int i = 1; i < N; i++) {
            value += ln(i + q) / doublePower(i + q, s) + ln(parameterHarmonicNumber(N, q, s));
        }
        return s / parameterHarmonicNumber(N, q, s) * value;
    }
}
