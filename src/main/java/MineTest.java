
import java.io.IOException;
import java.io.PrintStream;

import org.apache.commons.math3.exception.MathArithmeticException;
import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.util.FastMath;

public class MineTest {
    public static void main(String[] args) throws IOException {
        int a=114;
        int b=514;
        FastMath c=new FastMath();
        int d=c.min(a,b);
        System.out.println(d);
        int x=1919;
        int e=c.abs(x);
        System.out.println(x);
        double f=3;
        double g=c.sqrt(f);
        System.out.println(g);
    }
}