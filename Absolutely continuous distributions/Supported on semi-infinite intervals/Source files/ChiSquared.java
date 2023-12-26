public class ChiSquared extends MathLibrary{
    int k;

    ChiSquared(int k) {
        if (k <= 0) throw new IllegalArgumentException("k > 0");
        this.k = k;
    }

    double PDF(double x) {
        if (x < 0) throw new IllegalArgumentException("x >= 0");
        double part1 = 1 / (doublePower(2, k / 2.0) * gamma(k / 2.0));
        double part2 = doublePower(x, (k / 2.0) - 1) * euler(-1 * x / 2);
        return part1 * part2;
    }

    double CDF(double x) {
        if (x < 0) throw new IllegalArgumentException("x >= 0");
        return (1 / gamma(k / 2.0)) * incompleteLowerGamma(k / 2.0, x / 2);
    }

    int Mean() {
        return k;
    }

    double Median() {
        return k * power(1 - (2.0 / (9 * k)), 3);
    }

    double Mode() {
        return max(k - 2, 0);
    }

    int Variance() {
        return 2 * k;
    }

    double Skewness() {
        return root(8.0 / k, 2);
    }

    double EXkurtosis() {
        return 12.0 / k;
    }

    double Entropy() {
        double part1 = (k / 2.0) + ln(2 * gamma(k / 2.0));
        double part2 = (1 - (k / 2.0)) * digamma(k / 2.0);
        return part1 + part2;
    }

    double MGF(double t) {
        if (t < 0.5) {
            return doublePower(1 - (2 * t), -1 * k / 2.0);
        } else {
            return Double.NaN;
        }
    }

    double PGF(double t) {
        if (0 < t && t < root(euler(1), 2)) {
            return doublePower(1 - (2 * ln(t)), -1 * k / 2.0);
        } else {
            return Double.NaN;
        }
    }
}
