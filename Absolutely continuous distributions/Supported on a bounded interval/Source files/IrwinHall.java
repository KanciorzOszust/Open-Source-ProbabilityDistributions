import java.util.Random;

public class IrwinHall extends MathLibrary{
    int n;

    IrwinHall(int n) {
        if (n <= 0) throw new IllegalArgumentException("n > 0");
        this.n = n;
    }

    double PDF(double x) {
        double value = 1.0 / factorial(n - 1);
        for (int k = 0; k < (int) x; k++) {
            value += power(-1, k) * binomial(n, k) * power(x - k, n - 1);
        }
        return value;
    }

    double CDF(double x) {
        double value = 1.0 / factorial(n);
        for (int k = 0; k < (int) x; k++) {
            value += power(-1, k) * binomial(n, k) * power(x - k, n);
        }
        return value;
    }

    double Mean() {
        return n / 2.0;
    }

    double Median() {
        return Mean();
    }

    double Mode() {
        if (n == 1) {
            Random random = new Random();
            return random.nextDouble();
        } else {
            return Mean();
        }
    }

    double Variance() {
        return n / 12.0;
    }

    int Skewness() {
        return 0;
    }

    double EXkurtosis() {
        return -6.0 / (5 * n);
    }

    double MGF(double t) {
        return doublePower((euler(t) - 1) / t, t);
    }
}
