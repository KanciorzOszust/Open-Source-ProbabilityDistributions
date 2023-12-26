import java.util.Random;

public class Uniform extends MathLibrary{
    double a;
    double b;

    Uniform(double a, double b) {
        if (b > a) throw new IllegalArgumentException("a < b");
        this.a = a;
        this.b = b;
    }

    double PDF(double x) {
        if (x >= a && x <= b) {
            return 1 / (b - a);
        } else {
            return 0;
        }
    }

    double CDF(double x) {
        if (x < a) return 0;
        else if (x >= a && x <= b) {
            return (x - a) / (b - a);
        } else {
            return 1;
        }
    }

    double Mean() {
        return 0.5 * (a + b);
    }

    double Median() {
        return Mean();
    }

    double Mode() {
        Random random = new Random();
        return a + (b - a) * random.nextDouble();
    }

    double Variance() {
        return 1.0 / 12 * power(b - a, 2);
    }

    double MAD() {
        return 0.25 * (b - a);
    }

    int Skewness() {
        return 0;
    }

    double EXkurtosis() {
        return -6.0 / 5;
    }

    double Entropy() {
        return ln(b - a);
    }

    double MGF(double t) {
        if (t != 0) {
            return ((euler(t * a) - euler(t * b)) / (t * (b - a)));
        } else {
            return 1;
        }
    }
}
