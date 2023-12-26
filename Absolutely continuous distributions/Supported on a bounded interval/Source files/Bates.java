public class Bates extends MathLibrary{
    double a;
    double b;
    int n;

    Bates(double a, double b, int n) {
        if (b > a) throw new IllegalArgumentException("b > a");
        if (n < 1) throw new IllegalArgumentException(" n >= 1");
        this.a = a;
        this.b = b;
        this.n = n;
    }

    double PDF(double x) {
        double wartosc = (double) n / (2 * factorial(n - 1));
        for (int k = 0; k < n; k++) {
            wartosc += power(-1, k) * binomial(n, k) * power(n * x - k, n - 1) * sign(n * x - k);
        }
        return wartosc;
    }

    double Mean() {
        return 0.5 * (a + b);
    }

    double Variance() {
        return 1.0 / (12 * n) * power(b - a, 2);
    }

    int Skewness() {
        return 0;
    }

    double EXkurtosis() {
        return -6.0 / (5 * n);
    }
}
