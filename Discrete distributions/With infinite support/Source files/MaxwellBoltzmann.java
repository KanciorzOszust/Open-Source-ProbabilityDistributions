public class MaxwellBoltzmann extends MathLibrary {
    double a;

    MaxwellBoltzmann(double a) {
        if (a <= 0) throw new IllegalArgumentException("a > 0");
        this.a = a;
    }

    double PDF(double x) {
        double exponent = (- 1 * power(x, 2)) / (2 * power(a, 2));
        return root(2.0 / PI, 2) * (power(x, 2) / power(a, 3)) * euler(exponent);
    }

    double CDF(double x) {
        double exponent = (- 1 * power(x, 2)) / (2 * power(a, 2));
        double part1 = errorFunction((double) x / (root(2, 2) * a));
        double part2 = root(2.0 / PI, 2) * (x / a) * euler(exponent);
        return part1 - part2;
    }

    double Mean() {
        return 2 * a * root(2.0 / PI, 2);
    }

    double Mode() {
        return root(2, 2) * a;
    }

    double Variance() {
        return (double) (power(a, 2) * (3 * PI - 8)) / PI;
    }

    double Skewness() {
        double licznik = 2 * root(2, 2) * (16 - (5 * PI));
        double mianownik = (3 * PI - 8) * root(3 * PI - 8, 2);
        return licznik / mianownik;
    }

    double Entropy() {
        return ln(a * root(2 * PI, 2)) + EulerMascheroni - 0.5;
    }
}
