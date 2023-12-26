public class Nakagami extends MathLibrary{
    double m;
    double omega;

    Nakagami(double m, double omega) {
        if (m < 0.5) throw new IllegalArgumentException("m >= 0.5");
        if (omega <= 0) throw new IllegalArgumentException("omega > 0");
        this.m = m;
        this.omega = omega;
    }

    double PDF(double x) {
        if (x <= 0) throw new IllegalArgumentException("x > 0");
        double part1 = (2 * doublePower(m, m)) / (gamma(m) * doublePower(omega, m));
        double part2 = doublePower(x, 2 * m - 1) * euler(-1 * m / omega * power(x, 2));
        return part1 * part2;
    }

    double CDF(double x) {
        if (x <= 0) throw new IllegalArgumentException("x > 0");
        return incompleteLowerGamma(m, m / omega * power(x, 2)) / gamma(m);
    }

    double Mean() {
        return gamma(m + 0.5) / gamma(m) * root(omega / m, 2);
    }

    double Mode() {
        return (root(2, 2) / 2) * root(((2 * m - 1) * omega) / m, 2);
    }

    double Variance() {
        double part = (1 / m) * power(gamma(m + 0.5) / gamma(m), 2);
        return omega * (1 - part);
    }
}
