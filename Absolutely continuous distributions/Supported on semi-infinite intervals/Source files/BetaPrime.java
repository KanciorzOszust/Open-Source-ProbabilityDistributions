public class BetaPrime extends MathLibrary{
    double alpha;
    double beta;

    BetaPrime(double alpha, double beta) {
        if (alpha <= 0 || beta <= 0) throw new IllegalArgumentException("alpha, beta > 0");
        this.alpha = alpha;
        this.beta = beta;
    }

    double PDF(double x) {
        if (x < 0) throw new IllegalArgumentException("x >= 0");
        return (doublePower(x, alpha - 1) * doublePower(1 + x, -1 * alpha - beta)) / beta(alpha, beta);
    }

    double CDF(double x) {
        if (x < 0) throw new IllegalArgumentException("c >= 0");
        return incompleteBeta(alpha, beta, x / (1 + x));
    }

    double Mean() {
        if (beta > 1) return alpha / (beta - 1);
        else return Double.NaN;
    }

    double Mode() {
        if (alpha >= 1) return (alpha - 1) / (beta + 1);
        else return 0;
    }

    double Variance() {
        if (beta > 2) {
            double part1 = alpha * (alpha + beta - 1);
            double part2 = (beta - 2) * power(beta - 1, 2);
            return part1 / part2;
        } else {
            return Double.NaN;
        }
    }

    double Skewness() {
        if (beta > 3) {
            double part1 = (2 * (2 * alpha + beta - 1)) / (beta - 3);
            double part2 = root((beta - 2) / (alpha * (alpha + beta - 1)), 2);
            return part1 * part2;
        } else {
            return Double.NaN;
        }
    }
}
