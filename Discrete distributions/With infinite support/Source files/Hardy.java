public class Hardy extends MathLibrary{
    double p;
    double q;
    int m;

    Hardy(double p, double q, int m) {
        if (p <= 0 || p >= 1) throw new IllegalArgumentException("0 <= p <= 1");
        if (q <= 0 || q >= 1) throw new IllegalArgumentException("q <= q <= 1");
        if ((q + p) >= 1) throw new IllegalArgumentException("q + p < 1");
        if (m < 1) throw new IllegalArgumentException("m >= 1, n >= 1");
        this.p = p;
        this.q = q;
        this.m = m;
    }

    private double supportA(int j, int m) {
        int part1 = (int) (binomial(j - 1, 2 * j - m -1) * power(p, m - j + 1));
        int part2 = (int) (power(1 - p - q, 2 * j - m - 1));
        return part1 * part2;
    }

    private double supportB(int j, int m) {
        int part1 = (int) (binomial(j, 2 * j - m) * power(p, m - j));
        int part2 = (int) (power(1 - p - q, 2 * j - m));
        return part1 * part2;
    }

    double PMF(int n) {
        double wartosc = 0;
        int index = 0;
        if (m % 2 == 1) index = (m + 1) / 2;
        else index = m / 2;
        for (int i = index; i < m; i++) {
            wartosc += binomial(n - 1, n - i) * power(q, n - i) * (supportA(i, m) + supportB(i, m));
        }
        return wartosc;
    }

    double Mean() {
        double wartosc = 0;
        for (int i = 1; i < m; i++) {
            wartosc += ((m + 1 - i) * power(p, i - 1)) / power(q - 1, i);
        }
        return wartosc;
    }


    double MGF(double t) {
        double wartosc = 0;
        int index = 0;
        if (m % 2 == 1) index = (m + 1) / 2;
        else index = m / 2;
        for (int i = index; i < m; i++) {
            double licznik = (supportA(i, m) + supportB(i, m)) * euler(i * t);
            double mianownik = power(1 - (euler(t) * q), i);
            wartosc += licznik / mianownik;
        }
        return wartosc;
    }
}
