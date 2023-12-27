public class BetaBinominal extends MathLibrary {
    int n;
    double alpha;
    double beta;

    BetaBinominal(int n, double alpha, double beta) {
        if (n <= 0) throw new IllegalArgumentException("n > 0");
        if (alpha < 0) throw new IllegalArgumentException("alpha > 0");
        if (beta < 0) throw new IllegalArgumentException("beta > 0");
        this.n = n;
        this.alpha = alpha;
        this.beta = beta;
    }

    double PMF(int x) {
        if (x < 0 || x > n) throw new IllegalArgumentException("0 <= x <= n");
        return binomial(n, x) * (beta(x + alpha, n - x + alpha) / beta(alpha, beta));
    }

    double CDF(int x) {
        if (x < 0 || x > n) throw new IllegalArgumentException("0 <= x <= n");
        double value = 0;
        for (int i = 0; i < x; i++) {
            value += (gamma(alpha + i) * gamma(beta + n - i)) / (gamma(i + 1) * gamma(n - i + 1));
        }
        return (1 / beta(alpha, beta)) * (gamma(n + 1) / (gamma(alpha + beta + n))) * value;
    }

    double Mean() {
        return (n * alpha) / (alpha + beta);
    }

    double Variance() {
        double part1 = n * alpha * beta * (alpha + beta + n);
        double part2 = power(alpha + beta, 2) * (alpha + beta + 1);
        return part1 / part2;
    }

    double Skewness() {
        double part1 = (alpha + beta + (2 * n)) * (beta - alpha);
        double value1 = part1 / (alpha + beta + 2);
        double part2 = n * alpha * beta * (n + alpha + beta);
        double value2 = root((1 * alpha * beta) / part2, 2);
        return value1 * value2;
    }

    double MGF(double t) {
        return hypergeometric(-1 * n, alpha, alpha + beta, 1 - euler(t));
    }

    double PGF(double z) {
        return hypergeometric(-1 * n, alpha, alpha + beta, 1 - z);
    }
}
