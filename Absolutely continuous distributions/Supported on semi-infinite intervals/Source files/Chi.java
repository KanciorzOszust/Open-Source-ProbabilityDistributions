public class Chi extends MathLibrary{
    double k;

    Chi(double k) {
        if (k <= 0) throw new IllegalArgumentException("k > 0");
        this.k = k;
    }

    double PDF(double x) {
        if (x < 0) throw new IllegalArgumentException("x >= 0");
        double part1 = 1 / (doublePower(2, (k / 2) - 1) * gamma(k / 2));
        double part2 = doublePower(x, k - 1) * euler(-1 * power(x, 2) / 2);
        return part1 * part2;
    }

    double CDF(double x) {
        return regularizedGamma(k / 2, power(x, 2) / 2);
    }

    double Mean() {
        return root(2, 2) * gamma((k + 1) / 2) / gamma(k / 2);
    }

    double Median() {
        return root(k * power(1 - (2 / (9 * k)), 3), 2);
    }

    double Mode() {
        return root(k - 1, 2);
    }

    double Variance() {
        return k - root(Mean(), 2);
    }

    double Skewness() {
        double part1 = (Mean() / (Variance() * root(Variance(), 2)));
        double part2 = 1 - (2 * Variance());
        return part1 * part2;
    }

    double EXkurtosis() {
        return 2 / Variance() * (1 - (Mean() * root(Variance(), 2) * Skewness()) - Variance());
    }

    double Entropy() {
        double part1 = ln(gamma(k / 2));
        double part2 = k - ln(2) - ((k - 1) * polygamma(0, k / 2));
        return part1 + (0.5 * part2);
    }
}
