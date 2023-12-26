public class F extends MathLibrary{
    double d1;
    double d2;

    F(double d1, double d2) {
        if (d1 <= 0 || d2 <= 0) throw new IllegalArgumentException("d1,d2 > 0");
        this.d1 = d1;
        this.d2 = d2;
    }

    double PDF(double x) {
        if (x < 0) throw new IllegalArgumentException("x >= 0");
        double part1 = doublePower(d1 * x, d1) * doublePower(d2, d2);
        double part2 = doublePower(d1 * x + d2, d1 + d2);
        double part3 = x * beta(d1 / 2, d2 / 2);
        return root(part1 / part2, 2) / part3;
    }

    double CDF(double x) {
        if (x < 0) throw new IllegalArgumentException("x >= 0");
        return regularizedIncompleteBeta((d1 * x) / (d1 * x + d2), d1 / 2, d2 / 2);
    }

    double Mean() {
        if (d2 > 2) {
            return d2 / (d2 - 2);
        } else {
            return Double.NaN;
        }
    }

    double Mode() {
        if (d1 > d2) {
            return ((d1 - 2) / d1) * (d2 / (d2 + 2));
        } else {
            return Double.NaN;
        }
    }

    double variance() {
        if (d2 > 4) {
            double part1 = 2 * power(d2, 2) * (d1 + d2 - 2);
            double part2 = d1 * power(d2 - 2, 2) * (d2 - 4);
            return part1 / part2;
        } else {
            return Double.NaN;
        }
    }

    double Skewness() {
        if (d2 > 6) {
            double part1 = (2 * d1 + d2 - 2) * root(8 * (d2 - 4), 2);
            double part2 = (d2 - 6) * root(d1 * (d1 + d2 - 2), 2);
            return part1 / part2;
        } else {
            return Double.NaN;
        }
    }

    double EXkurtosis() {
        if (d2 > 8) {
            double part1 = (d1 * (5 * d2 - 22) * (d1 + d2 -2)) + ((d2 - 4) * power(d2 - 4, 2));
            double part2 = d1 * (d2 - 6) * (d2 - 8) * (d1 + d2 - 2);
            return 12 * (part1 / part2);
        } else {
            return Double.NaN;
        }
    }

    double Entropy() {
        return ln(gamma(d1 / d2)) + ln(gamma(d2 / 2)) - ln(gamma((d1 + d2) / 2))
        + (1 - (d1 / 2)) * digamma(1 + (d1 / 2)) - (1 + (d2 / 2)) * digamma(1 + (d2 / 2))
        + ((d1 + d2) / 2) * digamma((d1 + d2) / 2) + ln(d1 / d2);
    }
}
