public class Logistic extends MathLibrary{
    double mu;
    double s;

    Logistic(double mu, double s) {
        if (s <= 0) throw new IllegalArgumentException("s > 0");
        this.mu = mu;
        this.s = s;
    }

    double PDF(double x) {
        double part = euler(-1 * (x - mu) / s);
        return part / (s * power(1 + part, 2));
    }

    double CDF(double x) {
        return 1 / (1 + euler(-1 * (x - mu) / s));
    }

    double Quantile(double p) {
        return mu + (s * log(p / (1 - p), 10));
    }

    double Mean() {
        return mu;
    }

    double Median() {
        return mu;
    }

    double Mode() {
        return mu;
    }

    double Variance() {
        return power(s, 2) * power(PI, 2) / 3;
    }

    int Skewness() {
        return 0;
    }

    double EXkurtosis() {
        return 1.2;
    }

    double Entropy() {
        return ln(s) + 2;
    }

    double MGF(double t) {
        if (t > -1 / s && t < 1 / s) {
            return euler(mu * t) * beta(1 - (s * t), 1 + (s * t));
        } else {
            return Double.NaN;
        }
    }

    private double binaryEntropy(double x) {
        return -1 * x * ln(x) - ((1 - x) * ln(1 - x));
    }

    double ExpectedShortfall(double p) {
        return mu + ((s * binaryEntropy(p)) / (1 - p));
    }
}
