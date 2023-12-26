public class GeneralizedPareto extends MathLibrary{
    double mu;
    double sigma;
    double xi;

    GeneralizedPareto(double mu, double sigma, double xi) {
        if (sigma <= 0) throw new IllegalArgumentException("sigma > 0");
        this.mu = mu;
        this.sigma = sigma;
        this.xi = xi;
    }

    private double z(double x) {
        return (x - mu) / sigma;
    }

    double PDF(double x) {
        if (xi >= 0) {
            if (x < mu) throw new IllegalArgumentException("x >= mu");
        } else {
            if (x < mu || x > (mu - (sigma / xi))) throw new IllegalArgumentException("mu <= x <= mu - (sigma / xi)");
        }
        return (1 / sigma) * doublePower(1 + (xi * z(x)), -1 * 1 / (xi + 1));
    }

    double CDF(double x) {
        if (xi >= 0) {
            if (x < mu) throw new IllegalArgumentException("x >= mu");
        } else {
            if (x < mu || x > (mu - (sigma / xi))) throw new IllegalArgumentException("mu <= x <= mu - (sigma / xi)");
        }
        return 1 - doublePower(1 + (xi * z(x)), -1 / xi);
    }

    double Mean() {
        if (xi < 1) {
            return mu + (sigma / (1 - xi));
        } else {
            return Double.NaN;
        }
    }

    double Median() {
        return mu + (sigma * (doublePower(2, xi) - 1) / xi);
    }

    double Mode() {
        return mu;
    }

    double Variance() {
        if (xi < 0.5) {
            return power(sigma, 2) / (power(1 - xi, 2) * (1 - (2 * xi)));
        } else {
            return Double.NaN;
        }
    }

    double Skewness() {
        if (xi < 1.0 / 3) {
            double part = 2 * (1 + xi) * root(1 - (2 * xi), 2);
            return part / (1 - (3 * xi));
        } else {
            return Double.NaN;
        }
    }

    double EXkurtosis() {
        if (xi < 0.25) {
            double part1 = 3 * (1 - (2 * xi)) * (2 * power(xi, 2) + xi + 3);
            double part2 = (1 - (3 * xi)) * (1 - (4 * xi));
            return part1 / part2;
        } else {
            return Double.NaN;
        }
    }

    double Entropy() {
        return ln(sigma) + xi + 1;
    }

    private double MGFsupport(int j) {
        double value = 1;
        for (int k = 0; k < j; k++) {
            value *= 1 - (k * xi);
        }
        return value;
    }

    double MGF(double t) {
        double value = 0;
        for (int j = 0; j < 7; j++) {
            value += power(t * sigma, j) / MGFsupport(j);
        }
        return euler(t * mu) * value;
    }

    double ExpectedShortfall(double p) {
        double q = 1 - p;
        if (xi != 0) {
            double part = (doublePower(q, -1 * xi) / (1 - xi)) + ((doublePower(q, -1 * xi) - 1) / xi);
            return mu + (sigma * part);
        } else {
            return mu + (sigma * (1 - ln(q)));
        }
    }
}
