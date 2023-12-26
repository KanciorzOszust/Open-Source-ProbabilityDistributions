public class ExponentialLogarithmic extends MathLibrary{
    double p;
    double beta;

    ExponentialLogarithmic(double p, double beta) {
        if (p <= 0 || p >= 1) throw new IllegalArgumentException(" 0 < p < 1");
        if (beta <= 0) throw new IllegalArgumentException("beta > 0");
        this.p = p;
        this.beta = beta;
    }

    double PDF(double x) {
        if (x < 0) throw new IllegalArgumentException("x >= 0");
        double part1 = beta * (1 - p) * euler(-1 * beta * x);
        double part2 = 1 - (1 - p) * euler(-1 * beta * x);
        return (1 / (-1 * ln(p))) * (part1 / part2);
    }

    double CDF(double x) {
        if (x < 0) throw new IllegalArgumentException("x >= 0");
        return 1 - (ln(1 - (1 - p) * euler(-1 * beta * x)) / ln(p));
    }

    double Mean() {
        return -1 * polylogarithm(2, 1 - p) / (beta * ln(p));
    }

    double Median() {
        return ln(1 + root(p, 2)) / beta;
    }

    int Mode() {
        return 0;
    }

    double Variance() {
        return -1 * ((2 * polylogarithm(3, 1 - p)) / (power(beta, 2) * ln(p)));
    }

    double MGF(double t) {
        double part1 = -1 * (beta * (1 - p)) / ln(p * (beta - t));
        double part2 = hypergeometric(1, (beta - t) / beta, (2 * beta - t) / beta, 1 - p);
        return part1 * part2;
    }
}
