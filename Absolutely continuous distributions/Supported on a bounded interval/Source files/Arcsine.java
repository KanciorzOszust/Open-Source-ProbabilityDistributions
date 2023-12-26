import java.util.Random;

public class Arcsine extends MathLibrary{
    double PDF(double x) {
        return 1.0 / (PI * root(x * (1 - x), 2));
    }

    double CDF(double x) {
        return 2.0 / PI * arcSine(root(x, 2));
    }

    double Mean() {
        return 0.5;
    }

    double Median() {
        return 0.5;
    }

    int Mode() {
        Random random = new Random();
        return random.nextInt(1);
    }

    double Variance() {
        return 0.125;
    }

    int Skewness() {
        return 0;
    }

    double EXkurtosis() {
        return -1.5;
    }

    double Entropy() {
        return ln(PI / 4);
    }

    private double MGFpart(int k) {
        double wartosc = 1;
        for (int r = 0; r < k - 1; r++) {
            wartosc *= ((2 * r + 1) / (2 * r + 2));
        }
        return wartosc;
    }

    double MGF(double t) {
        double wartosc = 0;
        for (int k = 0; k < 7; k++) {
            wartosc += MGFpart(k) * power(t, k) * factorial(k);
        }
        return 1 + wartosc;
    }
}
