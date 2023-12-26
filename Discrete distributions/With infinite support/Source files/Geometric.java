public class Geometric extends MathLibrary {
    double p;

    Geometric(double p) {
        if (p <= 0 || p > 1) throw new IllegalArgumentException("0 < p <= 1");
        this.p = p;
    }

    double PMF(int k) {
        return power(1 - p, k - 1) * p;
    }

    double CDF(double x) {
        if (x < 1) return 0;
        else {
            return 1 - power(1 - p, (int) x);
        }
    }

    double Mean() {
        return 1 / p;
    }

    int Median() {
        return (int) ((-1 / log(1 - p, 2)) + 1);
    }

    int Mode() {
        return 1;
    }

    double Variance() {
        return (1 - p) / power(p, 2);
    }

    double Skewness() {
        return (2 - p) / root(1 - p, 2);
    }

    double ExKurtosis() {
        return 6 + (power(p, 2) / (1 - p));
    }

    double Entropy() {
        return (double) (-1 * (1 - p) * log(1 - p, 10) - (p * log(p, 10))) / p;
    }

    double MGF(double t) {
        if (t < (-1 * ln(1 - p))) {
            return (double) (p * euler(t)) / (1 - ((1 - p) * euler(t)));
        } else {
            return 0;
        }
    }

    double PGF(double z) {
        return (double) (p * z) / (1 - ((1 - p) * z));
    }
}
