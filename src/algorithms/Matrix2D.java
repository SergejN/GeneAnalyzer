/*
    File:
        Matrix2D.java
 *
    Revision:
        1.0.0.1
 *
    Description:
        Encapsulates the 2-dimensional matrix.
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package algorithms;

import java.awt.Dimension;
import java.awt.Point;


public class Matrix2D
{
    public static class EUnequalSizeException extends Exception
    {};

    public static class EIllegalSubscriptException extends Exception
    {};
    
    private int nRows = 0;
    private int nCols = 0;
    private float[][] nums = null;
    
    public Matrix2D(int nRows, int nColumns)
    {
        if(nRows<1 || nColumns<1)
            throw new IllegalArgumentException("Invalid matrix dimensions");
        this.nRows = nRows;
        this.nCols = nColumns;
        nums = new float[nRows][nColumns];
    }
    
    /**
     *  Returns the dimensions of the matrix.
     * 
     *  @return
     */
    public Dimension getDimensions()
    {
        return new Dimension(nRows, nCols);
    }
    
    /**
     *  Sets the value in the specified cell of the matrix.
     * 
     *  @param row
     *  @param column
     *  @param value
     */
    public void setValue(int row, int column, float value) throws EIllegalSubscriptException
    {
        if(row>=nRows || column>=nCols)
            throw new EIllegalSubscriptException();
        else
            nums[row][column]=value;
    }

    /**
     *  Returns the value in the specified cell of the matrix.
     *  If the coordinates are invalid, the method throws IllegalArgumentException.
     *
     *  Note:
     *  The first cell has the coordinates 0,0.
     *
     *  @param row
     *  @param column
     *  @return
     */
    public float getValue(int row, int column) throws EIllegalSubscriptException
    {
        if(row>=nRows || column>=nCols)
            throw new EIllegalSubscriptException();
        else
            return nums[row][column];
    }

    /**
     *  Returns the maximal value in the matrix.
     * 
     *  @return
     */
    public float getMaximum()
    {
        float max = nums[0][0];
        for(int i=0;i<nRows;i++)
        {
            for(int j=0;j<nCols;j++)
            {
                if(nums[i][j]>max)
                    max = nums[i][j];
            }
        }
        return max;
    }
    
    /**
     *  Returns the indices of the greatest matrix entry.
     * 
     *  @return
     */
    public Point getMaximumIndices()
    {
        Point p = new Point(0, 0);
        float max = nums[0][0];
        for(int i=0;i<nRows;i++)
        {
            for(int j=0;j<nCols;j++)
            {
                if(nums[i][j]>max)
                {
                    max = nums[i][j];
                    p.x = i;
                    p.y = j;
                }
            }
        }
        return p;
    }

    /**
     *  Returns the minimal value in the matrix.
     * 
     *  @return
     */
    public float getMinimum()
    {
        float min = nums[0][0];
        for(int i=0;i<nRows;i++)
        {
            for(int j=0;j<nCols;j++)
            {
                if(nums[i][j]<min)
                    min = nums[i][j];
            }
        }
        return min;
    }    

    /**
     *  Returns the indices of the smallest matrix entry.
     *
     *  @return
     */
    public Point getMinimumIndices()
    {
        Point p = new Point(0, 0);
        float min = nums[0][0];
        for(int i=0;i<nRows;i++)
        {
            for(int j=0;j<nCols;j++)
            {
                if(nums[i][j]<min)
                {
                    min = nums[i][j];
                    p.x = i;
                    p.y = j;
                }
            }
        }
        return p;
    }

    /**
     *  Returns the sum of the values in the matrix.
     * 
     *  @return
     */
    public float getSum()
    {
        float sum = 0.0f;
        for(int i=0;i<nRows;i++)
        {
            for(int j=0;j<nCols;j++)
            {
                sum += nums[i][j];
            }
        }
        return sum;
    }  
    
    /**
     *  Adds another matrix to the current and returns the resulting matrix.
     *  
     *  Note:
     *  Another matrix should have the same dimensions as the current one. 
     *  Otherwise the method throws an EUnequalSizeException. 
     * 
     *  @param m1
     *  @param m2
     *  @return
     */
    public Matrix2D add(Matrix2D matrix) throws EUnequalSizeException
    {
        if(nRows!=matrix.nRows || nCols!=matrix.nCols)
            throw new EUnequalSizeException();
        Matrix2D result = new Matrix2D(nRows, nCols);
        for(int i=0;i<nRows;i++)
        {
            for(int j=0;j<nCols;j++)
            {
                result.nums[i][j] = nums[i][j]+matrix.nums[i][j];
            }
        }
        return result;
    }
    
    /**
     *  Subtracts another matrix from the current and returns the resulting matrix.
     *  
     *  Note:
     *  Another matrix should have the same dimensions as the current one. 
     *  Otherwise the method throws an EUnequalSizeException. 
     * 
     *  @param m1
     *  @param m2
     *  @return
     */
    public Matrix2D subtract(Matrix2D matrix) throws EUnequalSizeException
    {
        if(nRows!=matrix.nRows || nCols!=matrix.nCols)
            throw new EUnequalSizeException();
        Matrix2D result = new Matrix2D(nRows, nCols);
        for(int i=0;i<nRows;i++)
        {
            for(int j=0;j<nCols;j++)
            {
                result.nums[i][j] = nums[i][j]-matrix.nums[i][j];
            }
        }
        return result;
    }

    /**
     *  Returns the submatrix of the current matrix.
     *
     *  @param firstRow
     *  @param lastRow
     *  @param firstCol
     *  @param lastCol
     *  @return
     */
    public Matrix2D submatrix(int firstRow, int lastRow, int firstCol, int lastCol) throws EIllegalSubscriptException
    {
        if(firstRow<0 || firstRow>=nRows || firstRow>lastRow || lastRow>=nRows ||
           firstCol<0 || firstCol>=nCols || firstCol>firstCol || lastCol>=nCols)
        {
            throw new EIllegalSubscriptException();
        }
        Matrix2D matrix = new Matrix2D(lastRow-firstRow+1, lastCol-firstCol+1);
        for(int i=firstRow;i<=lastRow;i++)
        {
            for(int j=firstCol;j<=lastCol;j++)
            {
                matrix.nums[i-firstRow][j-firstCol] = nums[i][j];
            }
        }
        return matrix;
    }

    /**
     *  Returns the transposed matrix.
     *
     *  @return
     */
    public Matrix2D transpose()
    {
        Matrix2D matrix = new Matrix2D(nCols, nRows);
        for(int i=0;i<=nRows;i++)
        {
            for(int j=0;j<=nCols;j++)
            {
                matrix.nums[j][i] = nums[i][j];
            }
        }
        return matrix;
    }

    @Override
    public String toString()
    {
        StringBuffer buf = new StringBuffer();
        for(int i=0;i<nRows;i++)
        {
            for(int j=0;j<nCols;j++)
            {
                buf.append(String.format("%f\t", nums[i][j]));
            }
            buf.append("\n");
        }
        return buf.toString();
    }
}
