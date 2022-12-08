package org.apache.commons.math3.util;


import static org.junit.Assert.*;

import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

public class FastMathTest {

    @BeforeClass
    public static void setUpBeforeClass() throws Exception {
    }

    @Test
    public void testSqrt() {

        double c = 4;
        double result;
        result = FastMath.sqrt(c);
        System.out.println(result);
        //fail("Not yet implemented");
    }

    @Test
    public void testAbsInt() {
        int c = -114;
        int result;
        result = FastMath.abs(c);
        System.out.println(result);
        //fail("Not yet implemented");
    }

    @Test
    public void testMinIntInt() {

        int c = 114,d = 514;
        int min1;
        min1 = FastMath.min(c,d);
        System.out.println(min1);
        //fail("Not yet implemented");
    }

}
