public class Benford extends MathLibrary{
    double PDF(int x) {
        if (x < 1) throw new IllegalArgumentException("x >= 1");
        return log(1 + (1 / x), 10);
    }

    double CDF(double x) {
        if (x < 1) throw new IllegalArgumentException("x >= 1");
        return log(1 + x, 10);
    }
}
