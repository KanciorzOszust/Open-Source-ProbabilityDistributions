public class Logarithmic extends MathLibrary {
    double p;

    Logarithmic(double p) {
        if (p <= 0 || p >= 1) throw new IllegalArgumentException("0 < p < 1");
        this.p = p;
    }

    double PMF(int k) {
        double part1 = -1.0 / ln(1 - p);
        double part2 = (double) power(p, k) / k;
        return part1 * part2;
    }

    double CDF(int k) {
        return 1 + ( (double) incompleteBeta(p, k + 1, 0) / ln(1 - p));
    }

    double Mean() {
        return (-1.0 / ln(1 - p)) * (p / (1 - p));
    }

    int Mode() {
        return 1;
    }

    double Variance() {
        double licznik = power(p, 2) + (p * ln(1 - p));
        double mianownik = power(1 - p, 2) * power(ln(1 - p), 2);
        return -1 * (licznik / mianownik);
    }

    double MGF(double t) {
        if (t < (-1 * ln(p))) {
            return (double) ln(1 - (p * euler(t))) / ln(1 - p);
        } else {
            return 0;
        }
    }

    double PGF(double z) {
        if (absolute(z) < 1.0 / p) {
            return (double) ln(1 - (p * z)) / ln(1 - p);
        } else {
            return 0;
        }
    }
}
