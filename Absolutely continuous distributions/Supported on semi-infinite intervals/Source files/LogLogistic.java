public class LogLogistic extends MathLibrary{
    double alpha;
    double beta;

    LogLogistic(double alpha, double beta) {
        if (alpha <= 0 || beta <= 0) throw new IllegalArgumentException("alpha, beta > 0");
        this.alpha = alpha;
        this.beta = beta;
    }

    double PDF(double x) {
        if (x < 0) throw new IllegalArgumentException("x >= 0");
        double part1 = (beta / alpha) * doublePower(x / alpha, beta - 1);
        double part2 = power(1 + doublePower(x / alpha, beta), 2);
        return part1 / part2;
    }

    double CDF(double x) {
        if (x < 0) throw new IllegalArgumentException("x >= 0");
        return 1 / (1 + doublePower(x / alpha, -1 * beta));
    }

    double Mean() {
        if (beta > 1) return (alpha * PI / beta) / sine(PI / beta);
        else return Double.NaN;
    }

    double Median() {
        return alpha;
    }

    double Mode() {
        if (beta > 1) {
            return alpha * doublePower((beta - 1) / (beta + 1), 1 / beta);
        } else {
            return 0;
        }
    }

    double Entropy() {
        return ln(alpha) - ln(beta) + 2;
    }

    double MGF(double t) {
        double value = 0;
        for (int n = 0; n < 7; n++) {
            value += power(alpha * t, n) / factorial(n) * beta(1 + (n / beta), 1 - (n / beta));
        }
        return value;
    }

    double ExpectedShortfall(double p) {
        double part1 = (PI / beta) * (1 / sine(PI / beta));
        double part2 = incompleteBeta(1 / beta + 1, 1 - (1 / beta), p);
        return (alpha / (1 - p)) * (part1 - part2);
    }
}
