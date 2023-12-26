public class BaldingNichols extends MathLibrary{
    double F;
    double p;
    double alpha;
    double beta;

    BaldingNichols(double F, double p) {
        if (F <= 0 || F >= 1 || p <= 0 || p >= 1) throw new IllegalArgumentException("0 < F, p < 1");
        this.F = F;
        this.p = p;
        this.alpha = ((1 - F) / F) * p;
        this.beta = ((1 - F) / F) * (1 - p);
    }

    double PDF(double x) {
        if (x <= 0 || x >= 1) throw new IllegalArgumentException("0 < x < 1");
        double part = doublePower(x, alpha - 1) * doublePower(1 - x, beta - 1);
        return part / beta(alpha, beta);
    }

    double CDF(double x) {
        if (x <= 0 || x >= 1) throw new IllegalArgumentException("0 < x < 1");
        return regularizedIncompleteBeta(alpha, beta, x);
    }

    double Mean() {
        return p;
    }

    double Median() {
        return 1 / regularizedIncompleteBeta(alpha, beta, 0.5);
    }

    double Mode() {
        return (F - ((1 - F) * p)) / (3 * F - 1);
    }

    double Variance() {
        return F * p * (1 - p);
    }

    double Skewness() {
        double part1 = 2 * F * (1 - (2 * p));
        double part2 = (1 + F) * root(F * (1 - p) * p, 2);
        return part1 / part2;
    }

    private double MGFsupport(double k) {
        double value = 1;
        for (int r = 0; r < k - 1; r++) {
            value *= (alpha + r) / (((1 - F) / F) + r);
        }
        return value;
    }

    double MGF(double t) {
        double value = 0;
        for (int k = 1; k < 7; k++) {
            value += MGFsupport(k) * power(t, k) / factorial(k);
        }
        return 1 + value;
    }
}
