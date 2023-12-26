public class Laplace extends MathLibrary{
    double mu;
    double b;

    Laplace(double mu, double b) {
        if (b <= 0) throw new IllegalArgumentException("b > 0");
        this.mu = mu;
        this.b = b;
    }

    double PDF(double x) {
        return (1 / (2 * b)) * euler(-1 * absolute(x - mu) / b);
    }

    double CDF(double x) {
        if (x <= mu) {
            return 0.5 * euler((x - mu) / b);
        } else {
            return 1 - (0.5 * euler(-1 * (x - mu) / b));
        }
    }

    double Quantile(double p) {
        if (p <= 0.5) {
            return mu + (b * ln(2 * p));
        } else {
            return mu - (b * ln(2 - (2 * p)));
        }
    }

    double Mean() {
        return mu;
    }

    double median() {
        return mu;
    }

    double Mode() {
        return mu;
    }

    double Variance() {
        return 2 * power(b, 2);
    }
    
    double MAD() {
        return b * ln(2);
    }

    int Skewness() {
        return 0;
    }

    int EXkurtosis() {
        return 3;
    }

    double Entropy() {
        return ln(2 * b * euler(1));
    }

    double MGF(double t) {
        if (absolute(t) < 1 / b) {
            return euler(mu * t) / (1 - (power(b, 2) * power(t, 2)));
        } else {
            return Double.NaN;
        }
    }

    double ExpectedShortfall(double p) {
        if (p < 0.5) {
            return mu + ((b * (p / (1 - p))) * (1 - ln(2 * p)));
        } else {
            return mu + (b * (1 - ln(2 * (1 - p))));
        }
    }
}
