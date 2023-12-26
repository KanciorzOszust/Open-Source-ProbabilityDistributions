public class Frechet extends MathLibrary{
    double a;
    double s = 1;
    double m = 0;

    Frechet(double a) {
        if (a <= 0) throw new IllegalArgumentException("a > 0");
        this.a = a;
    }

    Frechet(double a, double s, double m) {
        if (a <= 0 || s <= 0) throw new IllegalArgumentException("a,s > 0");
        this.a = a;
        this.s = s;
        this.m = m;
    }

    double PDF(double x) {
        if (x <= m) throw new IllegalArgumentException("x > m");
        double part1 = (a / s) * doublePower((x - m) / s, -1 - a);
        double part2 = euler(-1 * doublePower((x - m) / s, -1 * a));
        return part1 * part2;
    }

    double CDF(double x) {
        if (x <= m) throw new IllegalArgumentException("x > m");
        return euler(-1 * doublePower((x - m) / s, -1 * a));
    }

    double Mean() {
        if (a > 1) {
            return m + (s * gamma(1 - (1 / a)));
        } else {
            return Double.POSITIVE_INFINITY;
        }
    }

    double Median() {
        return m + (s / doublePower(ln(2), 1 / a));
    }

    double Mode() {
        return m + (s * doublePower(a / (1 + a), 1 / a));
    }

    double Variance() {
        if (a > 2) {
            return power(s, 2) * (gamma(1 - (2 / a)) - power(gamma(1 - (1 / a)), 2));
        } else {
            return Double.POSITIVE_INFINITY;
        }
    }

    double Skewness() {
        if (a > 3) {
            double part1 = gamma(1 - (4 / a)) - (3 * gamma(1 - (2 / a)) * gamma(1 / (1 / a))) + (2 * power(gamma(1 - (1 / a)), 3));
            double part2 = power(gamma(1 - (2 / a)) - power(gamma(1 - (1 / a)), 2), 3);
            return part1 / root(part2, 2);
        } else {
            return Double.POSITIVE_INFINITY;
        }
    }

    double EXkurtosis() {
        if (a > 4) {
            double part1 = gamma(1 - (4 / a)) - (3 * gamma(1 - (2 / a)) * gamma(1 - (1 / a))) + (3 * power(gamma(1 - (2 / a)), 2));
            double part2 = power(gamma(1 - (2 / a)) - power(gamma(1 - (1 / a)), 2), 2);
            return -6 + (part1 / part2);
        } else {
            return Double.POSITIVE_INFINITY;
        }
    }

    double Entropy() {
        return 1 + (EulerMascheroni / a) + EulerMascheroni + ln(s / a);
    }
}
