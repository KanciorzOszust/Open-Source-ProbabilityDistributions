public class GeneralizedNormal extends MathLibrary{
    double mu;
    double alpha;
    double beta;

    GeneralizedNormal(double mu, double alpha, double beta) {
        if (alpha <= 0 || beta <= 0) throw new IllegalArgumentException("alpha, beta > 0");
        this.mu = mu;
        this.alpha = alpha;
        this.beta = beta;
    }

    double PDF(double x) {
        double part1 = beta / (2 * alpha * gamma(1 / beta));
        double part2 = euler(-1 * doublePower(absolute(x - mu) / alpha, beta));
        return part1 * part2;
    }

    double CDF(double x) {
        double part1 = sign(x - mu) * (1 / (2 * gamma(1 / beta)));
        double part2 = incompleteLowerGamma(1 / beta, doublePower(absolute((x - mu) / alpha), beta));
        return part1 * part2 + 0.5;
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
        return power(alpha, 2) * gamma(3 / beta) / gamma(1 / beta);
    }

    int Skewness() {
        return 0;
    }

    double EXkurtosis() {
        return gamma(5 / beta) * gamma(1 / beta) / power(gamma(3 / beta), 2) - 3;
    }

    double Entropy() {
        return (1 / beta) - log(beta / (2 * alpha * gamma(1 / beta)), 10);
    }
}
