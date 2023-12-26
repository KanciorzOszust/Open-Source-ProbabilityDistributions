public class DiscreteUniform extends MathLibrary {
    int a;
    int b;
    int n;

    DiscreteUniform(int a, int b) {
        if (b < a) throw new IllegalArgumentException("b >= a");
        this.a = a;
        this.b = b;
        this.n = b - a + 1;
    }

    double PMF() {
        return 1.0 / n;
    }

    double CDF(int k) {
        return (double) (k - a + 1) / n;
    }

    double Mean() {
        return (a + b) / 2.0;
    }

    double Median() {
        return (a + b) / 2.0;
    }

    double Variance() {
        return (power(n, 2) - 1) / 12.0;
    }

    int Skewness() {
        return 0;
    }

    double EXkurtosis() {
        return -1 * (6 * (power(n, 2) + 1)) / (double) (5 * (power(n, 2) - 1));
    }

    double Entropy() {
        return ln(n);
    }

    double MGF(double t) {
        double licznik = euler(a * t) - euler((b + 1) * t);
        double mianownik = n * (1 - euler(t));
        return licznik / mianownik;
    }

    double PGF(double z) {
        double licznik = power(z, a) - power(z, b + 1);
        double mianownik = n * (1 - z);
        return licznik / mianownik;
    }
}


