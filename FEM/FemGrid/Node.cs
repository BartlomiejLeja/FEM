namespace FEM.FemGrid
{
    class Node
    {
        private double _x;
        private double _y;
        private bool _status;
        private double _temp;

        public Node(double x, double y, bool status, double temp)
        {
            _x = x;
            _y = y;
            _status = status;
            _temp = temp;
        }
            
        public double X { get=>_x; set=>_x=value; }
        public double Y { get=>_y; set=>_y=value; }
        public bool Status { get=>_status; set=>_status=value; }//true boundary conditions, false no boundary conditions
       public double Temp { get => _temp; set => _temp = 100; }
    }
}
