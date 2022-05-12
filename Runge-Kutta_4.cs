using System;
using System.Collections.Generic;
using System.Text;
using System.IO;
/*----------------------------------------------------------------------------------------------------------
|                                       Fourth-Order Runge-Kutta Solution for SIR Model                    |
|                                       Course:     AM6007                                                 |
|                                       Lecturer:   Dr Kieran Mulchrone                                    |
|                                                                                                          |
|                                       Name:       Mattis Voss                                            |
|                                       Student #:  121128764                                              |
|                                       Date:       15/11/2021                                             |
----------------------------------------------------------------------------------------------------------*/
namespace RungeKutta4
{
    class Program
    {
        static void Main(string[] args)
        {


            //---------------------------SIR MODEL--------------------------------------------------------

            // Set the parameters for the model and store them in a Vector
            double Rzero = 2.4;
            double gamma = 1 / 14.0;
            double beta = gamma * Rzero;
            Vector parameters = new Vector(new double[] { beta, gamma });

            // Create a FunctionVector of lambda functions (Here representing dS/dt, dI/dt, dR/dt).
            // First argument of FunctionVector constructor defines number of equations.
            FunctionVector f = new FunctionVector(3, parameters);
            f[0] = (x, y, theta) => -(theta[0] * y[0] * y[1]);
            f[1] = (x, y, theta) => (theta[0] * y[0] * y[1]) - (theta[1] * y[1]);
            f[2] = (x, y, theta) => theta[1] * y[1];

            // Create 4th order Runge-Kutta solver for the system, using the above FunctionVector as first argument.
            // Second argument is path to csv output file.
            RKSolver SIR = new RKSolver(f, "SIR.csv");

            // Set initial values for y (S = 0.99, I = 0.01, R = 0.00)
            Vector initValuesy = new Vector(new double[] { 0.99, 0.01, 0.00 });

            // Set initial and final values for x (i.e. time in this model)
            Vector valuesx = new Vector(new double[] { 0, 100 });

            // Set timestep
            double timeStep = 0.01;

            // Run solver
            SIR.Solve(initValuesy, valuesx, timeStep);


            // ---------------------------AM6005 CA2 MODEL-------------------------------------------

            // Set parameter values
            double mu = 1;
            double nu = 0;
            Vector para = new Vector(new double[] { mu, nu });

            // Define functions in terms of time t, variables y, and vector of parameters theta
            FunctionVector F = new FunctionVector(2, para);
            F[0] = (t, y, theta) => 2 * y[0] * (1 - y[0] / 2) - y[1] * (y[0] + theta[1]);
            F[1] = (t, y, theta) => y[1] * (theta[0] - y[1] * y[1]) - y[0] * (y[0] * y[1] - theta[1]);

            // New Solver for AM6005 
            RKSolver andy = new RKSolver(F, "AM6005.csv");

            // Set initial values for y
            Vector initValuesy = new Vector(new double[] { 0.01, 1 });

            // Set initial values for t
            Vector valuest = new Vector(new double[] { 0, 50 });

            // Set timestep
            double timeIncrement = 0.1;

            // Run solver
            andy.Solve(initValuesy, valuest, timeIncrement);




        }
    }



    /*=======================================================================================================================*/
    /*                                          RKSolver CLASS                                                               */
    /*=======================================================================================================================*/

    /// <summary>
    /// Solves a system of first-order ordinary differential equations. Takes in a FunctionVector of functions.
    /// Public methods:
    /// Solve():            Solves the system of equations for given initial values of y and x. Writes the solution
    ///                     to a csv file at each timestep.
    /// </summary>


    class RKSolver
    {
        /*-------------------------------------------------------------------------------------------------------------------*/
        /*                                          Private data                                                             */
        /*-------------------------------------------------------------------------------------------------------------------*/

        private FunctionVector f;
        private string path;

        /*-------------------------------------------------------------------------------------------------------------------*/
        /*                                          Constructors                                                             */
        /*-------------------------------------------------------------------------------------------------------------------*/

        public RKSolver(FunctionVector f, string path)
        {
            this.f = f;
            this.path = path;
        }

        // Solves the system of equations. First argument is a vector of initial values for y.
        // Second argument is the starting and ending values of x.
        // Third argument is the step size.

        public void Solve(Vector initValuesy, Vector valuesx, double stepSize)
        {

            if (valuesx[0] >= valuesx[1] || stepSize > (valuesx[1]-valuesx[0]) || valuesx[0] < 0)
            {
                Console.WriteLine("Please change time interval:");
                Console.WriteLine("xe should be larger than x0 >= 0. Step size should be smaller than time interval");
            }
            else
            {
                // Initialise all required variables for the Runge-Kutta algorithm
                double h = stepSize;
                double xj = valuesx[0];
                double xjp1;
                int numSteps = (int)((valuesx[1] - valuesx[0]) / stepSize);
                Vector yjp1 = new Vector(f.Elems.Length);
                Vector yj = new Vector(initValuesy);
                Vector k1 = new Vector(f.Elems.Length);
                Vector k2 = new Vector(f.Elems.Length);
                Vector k3 = new Vector(f.Elems.Length);
                Vector k4 = new Vector(f.Elems.Length);

                // Create a header line for a csv file
                string makeHeader(int length)
                {
                    string s = new string("x");

                    for (int i = 1; i <= length; i++)
                    {
                        s += ", y" + i.ToString();
                    }
                    return s;
                }
                
                // Create file to take output, "true" means existing information in file is not overwritten
                StreamWriter output = new StreamWriter(@path, true);

                // Check that file exists
                FileInfo finfo = new FileInfo(@path);
                try
                {
                    if (!finfo.Exists) throw new Exception("File does not exist");

                    // Create header line for csv file
                    string header = makeHeader(f.Elems.Length);
                    output.WriteLine(header);
                }
                catch (Exception)
                {
                    Console.WriteLine("{0}, cannot proceed");
                }

                // Fourth order Runge-Kutta iteration
                for (int j = 0; j <= numSteps; j++)
                {
                    output.WriteLine("{0}, {1}", xj, yj);

                    k1 = f.Eval(xj, yj);
                    k2 = f.Eval(xj + (h / 2), yj + (h / 2) * k1);
                    k3 = f.Eval(xj + (h / 2), yj + (h / 2) * k2);
                    k4 = f.Eval(xj + h, yj + h * k3);
                    yjp1 = yj + (1 / 6.0) * h * (k1 + 2 * k2 + 2 * k3 + k4);
                    xjp1 = xj + h;
                    xj = xjp1;
                    yj = yjp1;
                }
                output.Close();

                // Output values at final time step to the console
                Console.BackgroundColor = ConsoleColor.DarkBlue;
                Console.ForegroundColor = ConsoleColor.Yellow;
                Console.WriteLine("\n----Fourth-order Runge-Kutta solver-----");
                Console.WriteLine("\n----------------------------------------");
                Console.WriteLine("|   Final values for current solution  |");
                Console.WriteLine(string.Format("|  {0, -15}|{1, 18}  |", "x", (xj - h).ToString("0.000")));
                for (int i = 0; i < yj.Length; i++)
                    //Console.WriteLine("y{0}: {1}", i, yj[i]);
                    Console.WriteLine("|  y{0, -14}|{1, 18}  |", i+1, yj[i].ToString("0.00000"));
                Console.WriteLine("----------------------------------------\n");
                Console.ResetColor();
            }
        }
    }
    // Delegate for a method that takes a double and two vectors as input
    public delegate double func(double x, Vector y, Vector parameters);




    /*=======================================================================================================================*/
    /*                                          FunctionVector CLASS                                                         */
    /*=======================================================================================================================*/
    
    /// <summary>
    /// Creates a vector of functions that define a first-order ODE. The functions are assigned to a delegate of type 'func'.
    /// Has a public method 'Eval()', which evaluates each function in the FunctionVector and returns a vector of doubles.
    /// </summary>
    /// 
    public class FunctionVector
    {
        /*-------------------------------------------------------------------------------------------------------------------*/
        /*                                          Private data                                                             */
        /*-------------------------------------------------------------------------------------------------------------------*/
        // Array of function elements
        private func[] elems;
        // Vector of parameters of the equations
        private Vector theta;

        /*-------------------------------------------------------------------------------------------------------------------*/
        /*                                          Properties                                                               */
        /*-------------------------------------------------------------------------------------------------------------------*/
        public func[] Elems
        {
            get
            {
                return elems;
            }
            set
            {
                if (value.Length == this.Length) elems = value;
                else Console.WriteLine("Input array must have same length as existing vector");
            }
        }
        public int Length { get => elems.Length; }

        public Vector Parameters
        {
            get
            {
                return theta;
            }
            set
            {
                theta = value;
            }
        }
        /*-------------------------------------------------------------------------------------------------------------------*/
        /*                                          Constructors                                                             */
        /*-------------------------------------------------------------------------------------------------------------------*/
        
        // Constructor taking an array of functions and a vector of parameters
        public FunctionVector(func[] elems, Vector parameters)
        {
            this.elems = elems;
            this.theta = parameters;
        }
        // Constructor taking an integer defining number of functions as its argument
        public FunctionVector(int length)
        {
            func[] elems = new func[length];
            this.elems = elems;
        }

        // Constructor taking an integer defining number of functions and a vector of parameters
        public FunctionVector(int length, Vector parameters)
        {
            func[] elems = new func[length];
            this.elems = elems;
            this.theta = parameters;
        }

        /*-------------------------------------------------------------------------------------------------------------------*/
        /*                                          Indexing and public methods                                              */
        /*-------------------------------------------------------------------------------------------------------------------*/

        // Indexer for function array
        public func this[int index]
        {
            get { return this.elems[index]; }
            set { this.elems[index] = value; }
        }
        // Method to evaluate each function in the FunctionVector
        public Vector Eval(double x, Vector y)
        {
            Vector tmp = new Vector(this.elems.Length);
            for (int i = 0; i < this.elems.Length; i++)
                tmp[i] = this.Elems[i](x, y, this.theta);
            return tmp;
        }
    }


    /*=======================================================================================================================*/
    /*                                          Vector CLASS                                                                 */
    /*=======================================================================================================================*/

    /// <summary>
    /// Represents a vector. Length is fixed after initialistaion, but elements can be reassigned. 
    /// Allows use of [] notation for indexing.
    /// Default constructor asks user for input of vector length via console, creates a zero-Vector.
    /// 
    /// Overloaded operators:
    /// * Vector dot product (Vector * Vector)
    /// * Scaling of a vector by a double (double * Vector)
    /// + Binary and unary Vector addition
    /// - Binary and unary Vector subtraction
    /// 
    /// Public methods:
    /// ToString(): Returns vector as a row in square brackets.
    /// 
    /// Note that parts of this vector class were previously submitted by Mattis Voss 
    /// for assessment as part of the Perceptron assignment on 29/10/2021.
    /// 
    /// </summary>
    /// 
    public class Vector
    {
        /*-------------------------------------------------------------------------------------------------------------------*/
        /*                                          Private data                                                             */
        /*-------------------------------------------------------------------------------------------------------------------*/
        // Array of vector elements
        private double[] elems;
        /*-------------------------------------------------------------------------------------------------------------------*/
        /*                                          Properties                                                               */
        /*-------------------------------------------------------------------------------------------------------------------*/
        public double[] Elems
        {
            get
            {
                return elems;
            }
            set
            {
                if (value.Length == this.Length) elems = value;
                else Console.WriteLine("Input array must have same length as existing vector");
            }
        }
        public int Length { get => elems.Length; }

        /*-------------------------------------------------------------------------------------------------------------------*/
        /*                                          Constructors                                                             */
        /*-------------------------------------------------------------------------------------------------------------------*/

        // Default constructor, asks for length input from console
        public Vector()
        {
            int length;
            Console.WriteLine("Please enter length of vector required (positive integer):");
            string input = Console.ReadLine();
            bool success = int.TryParse(input, out length) && length >= 0;
            while (!success)
            {
                Console.WriteLine("Please try again, a positive integer is required:");
                input = Console.ReadLine();
                success = int.TryParse(input, out length) && length >= 0;
            }
            elems = new double[length];
        }
        // This constructor takes an array of doubles representing elements as its argument
        public Vector(double[] elems)
        {
            this.elems = elems;
        }
        // This constructor takes an integer defining length as its argument
        public Vector(int length)
        {
            bool success = length > 0 ? true : false;
            while (!success)
            {
                Console.WriteLine("Please enter a vector length greater than 0:");
                string input = Console.ReadLine();
                success = int.TryParse(input, out length) && length > 0;
            }

            double[] elems = new double[length];
            this.elems = elems;
        }
        // Copy constructor
        public Vector(Vector other)
        {
            this.elems = other.elems;
        }

        /*-------------------------------------------------------------------------------------------------------------------*/
        /*                                          Operator overloading                                                     */
        /*-------------------------------------------------------------------------------------------------------------------*/

        // Overload * operator for inner product
        public static double operator *(Vector left, Vector right)
        {
            double dotProd = 0;
            if (left.Length != right.Length)
                throw new ArgumentException("Inner product only possible for vectors of equal dimension");
            for (int i = 0; i < left.Length; i++)
            {
                dotProd += left.Elems[i] * right.Elems[i];
            }
            return dotProd;
        }
        // Overload + operator for binary addition
        public static Vector operator +(Vector left, Vector right)
        {
            Vector sum = new Vector(left.Elems.Length);
            if (left.Length != right.Length)
                throw new ArgumentException("Addition only possible for vectors of equal dimension");
            for (int i = 0; i < left.Length; i++)
            {
                sum.Elems[i] = left.Elems[i] + right.Elems[i];
            }
            return sum;
        }
        // Overload - operator for binary subtraction
        public static Vector operator -(Vector left, Vector right)
        {
            Vector diff = new Vector(left.Elems.Length);
            if (left.Length != right.Length)
                throw new ArgumentException("Subtraction only possible for vectors of equal dimension");
            for (int i = 0; i < left.Length; i++)
            {
                diff.Elems[i] = left.Elems[i] - right.Elems[i];
            }
            return diff;
        }

        // Overload + operator for unary addition
        public static Vector operator +(Vector right)
        {
            Vector tmp = new Vector(right.Elems.Length);
            for (int i = 0; i < right.Elems.Length; i++)
            {
                tmp.Elems[i] = right.Elems[i];
            }
            return tmp;
        }
        // Overload - operator for unary subtraction
        public static Vector operator -(Vector right)
        {
            int i = 0;
            Vector tmp = new Vector(right.Elems.Length);
            for (i = 0; i < right.Elems.Length; i++)
            {
                tmp.Elems[i] = -right.Elems[i];
            }
            return tmp;
        }
        // Overload * operator to scale a vector by a double
        public static Vector operator *(double x, Vector stretch)
        {
            Vector result = new Vector(stretch.Elems.Length);
            for (int i = 0; i < stretch.Length; i++)
            {
                result.Elems[i] = stretch.Elems[i] * x;
            }
            return result;
        }

        // Overload & operator to concatenate two vectors
        public static Vector operator &(Vector left, Vector right)
        {
            Vector tmp = new Vector(left.Length + right.Length);
            for (int i = 0; i < left.Length; i++)
            {
                tmp[i] = left[i];
            }
            for (int i = 0; i <= right.Length; i++)
            {
                tmp[left.Length + i + 1] = right[i];
            }
            return tmp;
        }
        /*-------------------------------------------------------------------------------------------------------------------*/
        /*                                          Indexing and overrides                                                   */
        /*-------------------------------------------------------------------------------------------------------------------*/

        // Indexer to allow use of [] notation
        public double this[int i]
        {
            get { return this.Elems[i]; }
            set { this.Elems[i] = value; }
        }
        // Override array string representation to return a row vector in square brackets
        public override string ToString()
        {
            string tmp = string.Join(", ", elems);
            return string.Format("{0}", tmp);
        }
        /*-------------------------------------------------------------------------------------------------------------------*/
        /*                                          Public methods                                                           */
        /*-------------------------------------------------------------------------------------------------------------------*/

        public double Norm()
        {
            double tmp = 0;
            for (int i = 0; i < 3; i++)
                tmp += elems[i] * elems[i];
            return Math.Sqrt(tmp);
        }
    }


}
