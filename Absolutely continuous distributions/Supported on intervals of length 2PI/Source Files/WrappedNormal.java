public class WrappedNormal extends MathLibrary{
    double mu;
    double sigma;

    WrappedNormal(double mu, double sigma) {
        if (sigma <= 0) throw new IllegalArgumentException("sigma > 0");
        this.mu = mu;
        this.sigma = sigma;
    }

    double PDF(double x) {
        if (x < -1 * PI || x > PI) throw new IllegalArgumentException("-PI <= x <= Pi");
        double value = 0;
        for (int i = -5; i < 5; i++) {
            double exponent = (-1 * power(x - mu + (2 * PI * i), 2)) / (2 * power(sigma, 2));
            value += euler(exponent);
        }
        return (1 / (sigma * root(2 * PI, 2))) * value;
    }

    double Mean() {
        if (mu > -1 * PI && mu < PI) return mu;
        else return Double.NaN;
    }

    double Median() {
        return Mean();
    }

    double Mode() {
        return mu;
    }

    double Variance() {
        return 1 - euler((-1 * power(sigma, 2)) / 2);
    }

    double Entropy() {
        double q = euler(-1 * power(sigma, 2));
        double value = 0;
        for (int k = 0; k < 7; k++) {
            value += (power(-1, k) / k) * (doublePower(1, (power(k, 2) + k) / 2) / (1 - power(q, k)));
        }
        double part = -1 * ln(eulerFunction(q) / (2 * PI));
        return part + (2 * value);
    }
}
