public class Dagum extends MathLibrary{
    double p;
    double a;
    double b;

    Dagum(double p, double a, double b) {
        if (a <= 0 || p <= 0 || b <= 0) throw new IllegalArgumentException("p,a,b > 0");
        this.p = p;
        this.a = a;
        this.b = b;
    }

    double PDF(double x) {
        if (x <= 0) throw new IllegalArgumentException("x > 0");
        double part = doublePower(doublePower(x / b, a) + 1, p + 1);
        return ((a * p) / x) * (doublePower(x / b, a * p) / part);
    }

    double CDF(double x) {
        if (x <= 0) throw new IllegalArgumentException(" x > 0");
        return doublePower(1 + doublePower(x / b, -1 * a), -1 * p);
    }

    double Mean() {
        if (a > 1) {
            return b * ((gamma(1 - (1 / a) * gamma(p + (1 / a)))) / gamma(p));
        } else {
            return Double.NaN;
        }
    }

    double Median() {
        return b * doublePower(-1 + doublePower(2, 1 / p), -1 / a);
    }

    double Mode() {
        return b * doublePower((a * p - 1) / (a + 1), 1 / a);
    }

    double Variance() {
        if (a > 2) {
            double part1 = (gamma(1 - (2 / a) * gamma(p + (2 / a)))) / gamma(p);
            double part2 = (gamma(1 - (1 / a)) * gamma(p + (1 / a))) / gamma(p);
            return power(b, 2) * (part1 - power(part2, 2));
        } else {
            return Double.NaN;
        }
    }
}
