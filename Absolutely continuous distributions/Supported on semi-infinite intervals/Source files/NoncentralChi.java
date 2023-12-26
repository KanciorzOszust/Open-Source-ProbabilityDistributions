public class NoncentralChi extends MathLibrary{
    double k;
    double gamma;

    NoncentralChi(double k, double gamma) {
        if (k <= 0 || gamma <= 0) throw new IllegalArgumentException("k, gamma > 0");
        this.k = k;
        this.gamma = gamma;
    }

    double PDF(double x) {
        if (x < 0) throw new IllegalArgumentException("x >= 0");
        double part1 = euler(-1 * (power(x, 2) + power(gamma, 2)) / 2) * doublePower(x, k) * gamma;
        double part2 = doublePower(gamma * x, k / 2);
        double part3 = modifiedBessel((k / 2) - 1, gamma * x);
        return part1 / part2 * part3;
    }

    double CDF(double x) {
        if (x < 0) throw new IllegalArgumentException("x >= 0");
        return 1 - marcumQ(k / 2, gamma, x);
    }

    double Mean() {
        return root(PI / 2, 2) * laguerre(0.5, (k / 2) - 1, -1 * root(gamma, 2) / 2);
    }

    double Variance() {
        return k + power(gamma, 2) - power(Mean(), 2);
    }
}
