public class GeneralizedExtremeValue extends MathLibrary{
    double mu;
    double sigma;
    double xi;

    GeneralizedExtremeValue(double mu, double sigma, double xi) {
        if (sigma <= 0) throw new IllegalArgumentException("sigma > 0");
        this.mu = mu;
        this.sigma = sigma;
        this.xi = xi;
    }

    private double support(double x) {
        if (xi != 0) {
            return doublePower(1 + (xi * ((x - mu) / sigma)), -1 / xi);
        } else {
            return euler(-1 * (x - mu) / sigma);
        }
    }

    double PDF(double x) {
        if (xi > 0) {
            if (x < ((mu - sigma) / xi)) throw new IllegalArgumentException("x >= (mu - sigma) / xi");
        } else if (xi < 0) {
            if (x > ((mu - sigma) / xi)) throw new IllegalArgumentException("x < (mu - sigma) / xi");
        }
        return (1 / sigma) * doublePower(support(x), xi + 1) * euler(-1 * support(x));
    }

    double CDF(double x) {
        if (xi > 0) {
            if (x < ((mu - sigma) / xi)) throw new IllegalArgumentException("x >= (mu - sigma) / xi");
        } else if (xi < 0) {
            if (x > ((mu - sigma) / xi)) throw new IllegalArgumentException("x < (mu - sigma) / xi");
        }
        return euler(-1 * support(x));
    }

    private double g(int k) {
        return gamma(1 - (k * xi));
    }

    double Mean() {
        if (xi < 0 && xi != 0) {
            return mu + ((sigma * (g(1) - 1)) / xi);
        } else if (xi == 0) {
            return mu + (sigma * EulerMascheroni);
        } else {
            return Double.POSITIVE_INFINITY;
        }
    }

    double Median() {
        if (xi != 0) {
            return mu + (sigma * ((doublePower(ln(2), -1 * xi) - 1) / xi));
        } else {
            return mu - (sigma * ln(ln(2)));
        }
    }

    double Variance() {
        if (xi < 0.5 && xi != 0) {
            return power(sigma, 2) * ((g(2) - power(g(1), 2)) / power(xi, 2));
        } else if (xi == 0) {
            return power(sigma, 2) * (power(PI, 2) / 6);
        } else {
            return Double.POSITIVE_INFINITY;
        }
    }

    double Skewness() {
        if (xi < 1 / 3) {
            double part1 = g(3) - (3 * g(2) * g(1)) + (2 * power(g(1), 3));
            double part2 = doublePower(g(2) * power(g(1), 2), 1.5);
            return sign(xi) * (part1 / part2);
        } else if (xi == 0) {
            return (12 * root(6, 2) * zeta(3)) / power(PI, 3);
        } else {
            return Double.NaN;
        }
    }

    double EXkurtosis() {
        if (xi < 0.25 && xi != 0) {
            double part = g(4) - (4 * g(3) * g(1)) + (6 * power(g(1), 2) * g(2)) - (3 * power(g(1), 4));
            return part / power(g(2) - power(g(1), 2), 2);
        } else if (xi == 0) {
            return 12.0 / 5;
        } else {
            return Double.NaN;
        }
    }

    double Entropy() {
        return ln(sigma) + (EulerMascheroni * xi) + EulerMascheroni + 1;
    }

    double ExpectedShortfall(double p) {
        if (xi != 0) {
            double part1 = mu + (sigma / (xi * (1 - p)));
            double part2 = incompleteLowerGamma(1 - xi, ln(1 / p)) - (1 - p);
            return part1 * part2;
        } else {
            double part1 = mu + (sigma / (1 - p));
            double part2 = -1 * p * ln(-1 * ln(p)) + logarithmicIntegral(p);
            return part1 * part2;
        }
    }
}
