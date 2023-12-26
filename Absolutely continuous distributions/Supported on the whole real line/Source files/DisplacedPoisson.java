public class DisplacedPoisson extends MathLibrary{
    double lambda;
    double r;

    DisplacedPoisson(double lambda, double r) {
        if (lambda <= 0) throw new IllegalArgumentException("lambda > 0");
        this.lambda = lambda;
        this.r = r;
    }

    double Mean() {
        return lambda - r;
    }

    int Mode() {
        if (lambda >= (r + 1)) {
            return (int) (lambda - r);
        } else {
            return 0;
        }
    }

    double Variance() {
        return lambda;
    }

    private double MGFspport(double x, double y) {
        double wartosc = 0;
        for (int i = (int) x; i < 7; i++) {
            wartosc += (euler(-1 * y) * power(y, i)) / factorial(i);
        }
        return wartosc;
    }

    double MGF(double t, int k) {
        double wartosc = euler((lambda * euler(t - 1)) - (t * r));
        if (r >= 0) {
            double part2 = MGFspport(r + k, lambda * euler(t)) / MGFspport(r + k, lambda);
            wartosc *= part2;
        }
        return wartosc;
    }
}
