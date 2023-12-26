public class Levy extends MathLibrary{
    double mu;
    double c;

    Levy(double mu, double c) {
        if (c <= 0) throw new IllegalArgumentException("c > 0");
        this.mu = mu;
        this.c = c;
    }

    double PDF(double x) {
        if (x <= mu) throw new IllegalArgumentException("x > mu");
        double part1 = root(x / (2 * PI), 2);
        double part2 = euler(-1 * c / (2 * (x - mu))) / ((x - mu) * root(x - mu, 2));
        return part1 * part2;
    }

    double CDF(double x) {
        if (x <= mu) throw new IllegalArgumentException("x > mu");
        return ERFC(root(c / (2 * (x - mu)), 2));
    }

    double Quantile(double p) {
        return mu + (c / (2 * power(INVERFC(p), 2)));
    }

    double Mean() {
        return Double.POSITIVE_INFINITY;
    }

    double Median() {
        return mu + ((c / 2) * power(INVERFC(0.5), 2));
    }

    double Variance() {
        return Double.POSITIVE_INFINITY;
    }

    double Entropy() {
        return (1 + (3 * EulerMascheroni) + ln(16 * PI * power(c, 2))) / 2;
    }
}
