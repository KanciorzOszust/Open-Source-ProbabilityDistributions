public class WignerSemicircle extends MathLibrary{
    double R;
    double squareR;

    WignerSemicircle(double R) {
        if (R <= 0) throw new IllegalArgumentException("R > 0");
        this.R = R;
        this.squareR = power(R, 2);
    }

    double PDF(double x) {
        if (x < -1 * R || x > R) throw new IllegalArgumentException("-R <= x <= R");
        return (2 / (PI * squareR)) * root(squareR - power(x, 2), 2);
    }

    double CDF(double x) {
        if (x < -1 * R || x > R) throw new IllegalArgumentException("-R <= x <= R");
        double part1 = ((x * root(squareR - power(x, 2), 2)) / (PI * squareR));
        double part2 = arcSine(x / R) / PI;
        return 0.5 + part1 + part2;
    }

    int Mean() {
        return 0;
    }

    int Median() {
        return 0;
    }

    int Mode() {
        return 0;
    }

    double Variance() {
        return squareR / 4;
    }

    int Skewness() {
        return 0;
    }

    int EXkurtosis() {
        return -1;
    }

    double Entropy() {
        return ln(PI * R) - 0.5;
    }

    double MGF(double t) {
        return 2 * (modifiedBessel(1, R * t) / (R * t));
    }
}
