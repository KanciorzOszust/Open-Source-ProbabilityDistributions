public class Hermite extends MathLibrary{
    double a1;
    double a2;

    Hermite(double a1, double a2) {
        if (a1 < 0 || a2 < 0) throw new IllegalArgumentException("a1 >= 0, a2 >= 0");
        this.a1 = a1;
        this.a2 = a2;
    }

    double PMF(int x) {
        double wartosc = euler(-1 * (a1 + a2));
        for (int i = 0; i < (int) x / 2.0; i++) {
            double licznik = power(a1, x - 2 * i) * power(a2, i);
            double mianownik = factorial(x - 2 * i) * factorial(i);
            wartosc += licznik / mianownik;
        }
        return wartosc;
    }

    double CDF(int x) {
        double wartosc = euler(-1 * a1 + a2);
        for (int i = 0; i < x; i++) {
            for (int j = 0; j < (int) i / 2.0; j++) {
                double licznik = power(a1, i - 2 * j) * power(a2, j);
                double mianownik = factorial(i - 2 * j) * factorial(j);
                wartosc += licznik / mianownik;
            }
        }
        return wartosc;
    }

    double Mean() {
        return a1 + (2 * a2);
    }

    double Variance() {
        return a1 + (4 * a2);
    }

    double Skewness() {
        double licznik = a1 + (8 * a2);
        double mianownik = (a1 + (4* a2)) * root(a1 + (4 * a2), 2);
        return licznik / mianownik;
    }

    double EXkurtosis() {
        return (a1 + (16 * a2)) / power(a1 + (4 * a2), 2);
    }

    double MGF(double t) {
        return euler(a1 * (euler(t) - 1) + a2 * (euler(2 * t) - 1));
    }

    double PGF(double s) {
        return euler(a1 * (s - 1) + a2 * (power(s, 2) - 1));
    }
}
