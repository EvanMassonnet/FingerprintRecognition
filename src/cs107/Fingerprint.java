package cs107;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Provides tools to compare fingerprint.
 */
public class Fingerprint {

  /**
   * The number of pixels to consider in each direction when doing the linear
   * regression to compute the orientation.
   */
  public static final int ORIENTATION_DISTANCE = 16;

  /**
   * The maximum distance between two minutiae to be considered matching.
   */
  public static final int DISTANCE_THRESHOLD = 5;

  /**
   * The number of matching minutiae needed for two fingerprints to be considered
   * identical.
   */
  public static final int FOUND_THRESHOLD = 20;

  /**
   * The distance between two angle to be considered identical.
   */
  public static final int ORIENTATION_THRESHOLD = 20;

  /**
   * The offset in each direction for the rotation to test when doing the
   * matching.
   */
  public static final int MATCH_ANGLE_OFFSET = 2;

  /**
   * Returns an array containing the value of the 8 neighbours of the pixel at
   * coordinates <code>(row, col)</code>.
   * <p>
   * The pixels are returned such that their indices corresponds to the following
   * diagram:<br>
   * ------------- <br>
   * | 7 | 0 | 1 | <br>
   * ------------- <br>
   * | 6 | _ | 2 | <br>
   * ------------- <br>
   * | 5 | 4 | 3 | <br>
   * ------------- <br>
   * <p>
   * If a neighbours is out of bounds of the image, it is considered white.
   * <p>
   * If the <code>row</code> or the <code>col</code> is out of bounds of the
   * image, the returned value should be <code>null</code>.
   *
   * @param image array containing each pixel's boolean value.
   * @param row   the row of the pixel of interest, must be between
   *              <code>0</code>(included) and
   *              <code>image.length</code>(excluded).
   * @param col   the column of the pixel of interest, must be between
   *              <code>0</code>(included) and
   *              <code>image[row].length</code>(excluded).
   * @return An array containing each neighbours' value.
   */
  public static boolean[] getNeighbours(boolean[][] image, int row, int col) {
	  assert (image != null);

	  final int width = image[0].length;		//col
	  final int height = image.length;			//row
	  
	  if (row < 0 || row >= height || col < 0 || col >= width) {
		  return null;
	  }

	  boolean[] neighbours = new boolean[8];	//tableau de 8 elements qui sera renvoyer

	  //on remplie manuellement les 8 cases du tableau en verifiant a chaque fois que l'emplacement du tableau exist
	  // Si il n'existe pas c'est un rebord du tableau et on renvoie false par defaud
	  neighbours[0] = (row > 0)?							image[row-1][col] : false;
	  neighbours[1] = (row > 0 && col < width-1)? 			image[row-1][col+1] : false;
	  neighbours[2] = (col < width-1)? 						image[row][col+1] : false;
	  neighbours[3] = (row < height-1 && col < width-1)?	image[row+1][col+1] : false;
	  neighbours[4] = (row < height-1)? 					image[row+1][col] : false;
	  neighbours[5] = (row < height-1 && col > 0)? 			image[row+1][col-1] : false;
	  neighbours[6] = (col > 0)? 							image[row][col-1] : false;
	  neighbours[7] = (row > 0 && col > 0)? 				image[row-1][col-1] : false;
	  
	  return neighbours;
  }

  /**
   * Computes the number of black (<code>true</code>) pixels among the neighbours
   * of a pixel.
   *
   * @param neighbours array containing each pixel value. The array must respect
   *                   the convention described in
   *                   {@link #getNeighbours(boolean[][], int, int)}.
   * @return the number of black neighbours.
   */
  public static int blackNeighbours(boolean[] neighbours) {

	  assert(neighbours.length == 8);

	  int nbNeighbours = 0;		//variable qui contient le nombre de voisin et qui sera retourner
	  
	  for(int i = 0; i < 8; ++i) {
		  nbNeighbours += (neighbours[i]) ? 1 : 0;
	  }
	  
	  return nbNeighbours;
  }
  
  /**
   * Computes the number of white to black transitions among the neighbours of
   * pixel.
   *
   * @param neighbours array containing each pixel value. The array must respect
   *                   the convention described in
   *                   {@link #getNeighbours(boolean[][], int, int)}.
   * @return the number of white to black transitions.
   */
  public static int transitions(boolean[] neighbours) {

	  assert(neighbours.length == 8);
	  
	  int nbTransition = 0;		//variable qui contient le nombre de transition et qui sera renvoyer
	  for(int i =0 ; i < 7; ++i ) {
		  if(!neighbours[i] && neighbours[i+1]) {
			  ++nbTransition;
		  }
	  }
	  //on fini par comparer la transition du dernier element et du premier
	  if(!neighbours[7] && neighbours[0]) {
		  ++nbTransition;
	  }

	  return nbTransition;
  }

  /**
   * Returns <code>true</code> if the images are identical and false otherwise.
   *
   * @param image1 array containing each pixel's boolean value.
   * @param image2 array containing each pixel's boolean value.
   * @return <code>True</code> if they are identical, <code>false</code>
   *         otherwise.
   */
  public static boolean identical(boolean[][] image1, boolean[][] image2) {
	  if(image1 == null || image2 ==null) {		//si l'une des deux image est null elles ne peuvent pas etre egale
		  return false;
	  }
	  
	  if(image1.length != image2.length || image1[0].length != image2[0].length) { //si leur taille sont differente elle ne peuvent pas etre egale
		  return false;
	  }

	  final int width = image1[0].length;
	  final int height = image1.length;
	  int row =0;
	  int col =0;
	  
	  
	  while(row < height) {			//on parcours les iamge suivant leur lignes
		  col=0;					//chaque fois qu'on change de colonne on reinitialise la variable col a 0
		  while(col < width) {			//on parcours les iamge suivant leur colonnes
			  if(image1[row][col] != image2[row][col]) {
				  return false;									//si un seul element est different on renvoie imediatement false et la fonction se termine
			  }
			  ++col;
		  }
		  ++row;
	  }
	  return true;				//finalement si tout les element sont identique on renvoie true
  }


	/**
	 * Return a copy of a 2D array
     *
	 * @param array array of 2 dimensions
	 * @return a copy of the original 2D array
	 */
  public static boolean[][] copy2DArray(boolean[][] array) {
	  assert (array != null);
	  
	  final int width = array[0].length;
	  final int height = array.length;
	  
	  boolean[][] newImage = new boolean[height][width];        //on initialise un tableau vide de meme taille que celui placé en entrer
	  
	  for(int i =0; i < height; ++i) {
		  System.arraycopy(array[i], 0, newImage[i], 0, width);     //on utilise la fonction java qui permet de copier un tableau de dimension 1 sur chaques colonnes
	  }
	  return newImage;
  }


  /**
   * Internal method used by {@link #thin(boolean[][])}.
   *
   * @param image array containing each pixel's boolean value.
   * @param step  the step to apply, Step 0 or Step 1.
   * @return A new array containing each pixel's value after the step.
   */
  public static boolean[][] thinningStep(boolean[][] image, int step) {
	  
	  assert (image != null);
      assert (step == 0 || step ==1);
	  
	  final int width = image[0].length;
	  final int height = image.length;

	  boolean[][] newImage = copy2DArray(image);        //on creéé une copie du tableau placé en entrer pour eviter de le modifier
	  
	  for(int row =0; row < height; ++row) {        //on parcours chaque valeur du tableau
		  for(int col=0; col < width; ++col) {

			  boolean[] neighboursList = getNeighbours(image, row, col);

			  //chaque condition est verifié sur l'image non modifier

			  if(	image[row][col] &&                                  //si le pixel est noir
					neighboursList.length != 0  &&                      //si le pixel a des voisins
					blackNeighbours(neighboursList) >=2 && blackNeighbours(neighboursList) <=6 &&       //si il a entre 2 et 6 voisins noir
					transitions(neighboursList) == 1) {                                                 //et une unique transition

				  //le reste des conditions depend de l'etape
				  if( step == 0 &&
						  (!neighboursList[0] || !neighboursList[2] || !neighboursList[4]) &&
							(!neighboursList[2] || !neighboursList[4] || !neighboursList[6])) {
					  
					  newImage[row][col] = false;                           //on supprime le pixel noir (on le remplace par un blanc)
					  
				  }else if (step == 1 &&
						  (!neighboursList[0] || !neighboursList[2] || !neighboursList[6]) &&
							(!neighboursList[0] || !neighboursList[4] || !neighboursList[6])) {

					  newImage[row][col] = false;
					  
				  }	  
			   }  
			  }
		  }
	  
	  return newImage;      //on renvoie la copie modifier
	  }


  public static boolean[][] thinning(boolean[][] image){
	  boolean[][] newImage = thinningStep(image, 0);
	  newImage = thinningStep(newImage, 1);
	  return newImage;
  }
  /**
   * Compute the skeleton of a boolean image.
   *
   * @param image array containing each pixel's boolean value.
   * @return array containing the boolean value of each pixel of the image after
   *         applying the thinning algorithm.
   */
  public static boolean[][] thin(boolean[][] image) {
	  
	  assert (image != null);

	  boolean[][] previousImage = copy2DArray(image);
	  boolean[][] nextImage = thinningStep(image, 0);
	  nextImage = thinningStep(nextImage, 1);
	  
	  while(!identical(previousImage, nextImage)){
		  previousImage = nextImage;
		  nextImage = thinningStep(nextImage, 0);       //on applique les deux transformations
		  nextImage = thinningStep(nextImage, 1);
	  }
	  return nextImage;
  }


  /**
   * Computes all pixels that are connected to the pixel at coordinate
   * <code>(row, col)</code> and within the given distance of the pixel.
   *
   * @param image    array containing each pixel's boolean value.
   * @param row      the first coordinate of the pixel of interest.
   * @param col      the second coordinate of the pixel of interest.
   * @param distance the maximum distance at which a pixel is considered.
   * @return An array where <code>true</code> means that the pixel is within
   *         <code>distance</code> and connected to the pixel at
   *         <code>(row, col)</code>.
   */
  public static boolean[][] connectedPixels(boolean[][] image, int row, int col, int distance) {

	  assert (image != null);
	  assert (distance >= 0);

	  final int height = image.length;
	  final int width = image[0].length;

	  boolean[][] connectedPixelsImage = new boolean[height][width];

	  if ( !image[row][col] || row < 0 || row >= image.length || col < 0 || col >= image[0].length){
		  return connectedPixelsImage;
	  }

	  ArrayList<Integer> pixelsToDoRow = new ArrayList<Integer>();
	  ArrayList<Integer> pixelsToDoCol = new ArrayList<Integer>();

	  pixelsToDoRow.add(row);
	  pixelsToDoCol.add(col);

	  int actualRow;
	  int actualCol;

	  while(!pixelsToDoRow.isEmpty()){

		  actualRow = pixelsToDoRow.get(0);
		  actualCol = pixelsToDoCol.get(0);

		  pixelsToDoRow.remove(0);
		  pixelsToDoCol.remove(0);

		  for(int i = actualRow - 1; i < actualRow +2; ++i){
			  for(int j = actualCol - 1; j < actualCol +2; ++j){
				  if( i >= 0 && i < height && j >= 0 && j < width){
					  if(image[i][j] && !connectedPixelsImage[i][j] && !(i == actualRow && j == actualCol)){
						  if(Math.abs(row-i) <= distance && Math.abs(col-j) <= distance){
							  pixelsToDoRow.add(i);
							  pixelsToDoCol.add(j);
						  }
					  }
				  }
			  }
		  }
		  connectedPixelsImage[actualRow][actualCol] = true;
	  }

	  return connectedPixelsImage;
  }

  /**
   * Computes the slope of a minutia using linear regression.
   *
   * @param connectedPixels the result of
   *                        {@link #connectedPixels(boolean[][], int, int, int)}.
   * @param row             the row of the minutia.
   * @param col             the col of the minutia.
   * @return the slope.
   */
  public static double computeSlope(boolean[][] connectedPixels, int row, int col) {

	  assert (connectedPixels != null);
	  assert (connectedPixels[row][col]);

	  //col = x
	  //row = y
	  final int height = connectedPixels.length;
	  final int width = connectedPixels[0].length;

	  double sumXSquare = 0.0;
	  double sumYSquare = 0.0;
	  double sumProduct = 0.0;

	  for(int actualRow =0; actualRow < height; ++actualRow) {
		  for(int actualCol=0; actualCol < width; ++actualCol) {
			  if(connectedPixels[actualRow][actualCol]){
				sumXSquare += Math.pow((actualCol - col), 2);
				sumYSquare += Math.pow((row - actualRow), 2);
				sumProduct += (actualCol - col) * (row - actualRow);
			  }
		  }
	  }
	  if(sumXSquare == 0){
		  return Double.POSITIVE_INFINITY;
	  }

	  if(sumXSquare >= sumYSquare){
		  return sumProduct / sumXSquare;
	  }

	  return sumYSquare / sumProduct;
  }

  /**
   * Computes the orientation of a minutia in radians.
   * 
   * @param connectedPixels the result of
   *                        {@link #connectedPixels(boolean[][], int, int, int)}.
   * @param row             the row of the minutia.
   * @param col             the col of the minutia.
   * @param slope           the slope as returned by
   *                        {@link #computeSlope(boolean[][], int, int)}.
   * @return the orientation of the minutia in radians.
   */
  public static double computeAngle(boolean[][] connectedPixels, int row, int col, double slope) {

	  assert (connectedPixels != null);
	  assert (connectedPixels[row][col]);

	  final int height = connectedPixels.length;
	  final int width = connectedPixels[0].length;

	  int top = 0;
	  int bottom = 0;

      if(slope == Double.POSITIVE_INFINITY){
          for(int actualRow =0; actualRow < height; ++actualRow) {
                  if(actualRow != row && connectedPixels[actualRow][col]){
                      if((row - actualRow) > 0){
                          ++top;
                      }else{
                          ++bottom;
                      }
                  }
          }
		  return (top >= bottom) ? Math.PI/2 : -Math.PI/2;
      }

	  double perpSlope = -1/slope;

	  for(int actualRow =0; actualRow < height; ++actualRow) {
		  for(int actualCol=0; actualCol < width; ++actualCol) {
			  if((actualRow != row || actualCol != col) && connectedPixels[actualRow][actualCol]){
				if((row - actualRow) >= perpSlope * (actualCol - col)){
					++top;
				}else{
					++bottom;
				}
			  }
		  }
	  }

	  if((Math.atan(slope) >= 0 && bottom > top) || (Math.atan(slope) < 0 && bottom < top)){
		  return Math.atan(slope) + Math.PI;
	  }else{
		  return Math.atan(slope);
	  }
  }

  /**
   * Computes the orientation of the minutia that the coordinate <code>(row,
   * col)</code>.
   *
   * @param image    array containing each pixel's boolean value.
   * @param row      the first coordinate of the pixel of interest.
   * @param col      the second coordinate of the pixel of interest.
   * @param distance the distance to be considered in each direction to compute
   *                 the orientation.
   * @return The orientation in degrees.
   */
  public static int computeOrientation(boolean[][] image, int row, int col, int distance) {


      boolean[][] returnImage = connectedPixels(image, row, col, distance);
      double slope = computeSlope(returnImage,row, col);
      double angle = computeAngle(returnImage, row, col, slope);

      int angleDegrees = (int)Math.round(Math.toDegrees(angle));
      if(angleDegrees < 0){
          return angleDegrees+360;
      }
      return  angleDegrees;
  }


  public static List<int[]> copyMinutiaeList(List<int[]> minutiae) {

		List<int[]> newList = new ArrayList<int[]>();

		for(int m = 0; m < minutiae.size(); ++m){

			int[]minutia = {minutiae.get(m)[0], minutiae.get(m)[1], minutiae.get(m)[2]};
			newList.add(minutia);
		}
		return newList;
	}

  /**
   * Extracts the minutiae from a thinned image.
   *
   * @param image array containing each pixel's boolean value.
   * @return The list of all minutiae. A minutia is represented by an array where
   *         the first element is the row, the second is column, and the third is
   *         the angle in degrees.
   * @see #thin(boolean[][])
   */
  public static List<int[]> extract(boolean[][] image) {

	  assert (image != null);

	  List<int[]> minutia = new ArrayList<int[]>();

	  final int height = image.length;
	  final int width = image[0].length;

	  for(int row =1; row < height-1; ++row) {
		  for(int col=1; col < width-1; ++col) {
			  if(image[row][col]){
				  if(transitions(getNeighbours(image,row,col)) == 1 || transitions(getNeighbours(image,row,col)) ==3){
					  int[] minutiafound = {row, col,computeOrientation(image, row, col, ORIENTATION_DISTANCE)};
					  minutia.add(minutiafound);
				  }
			  }
		  }
	  }

	  return minutia;
  }

  /**
   * Applies the specified rotation to the minutia.
   *
   * @param minutia   the original minutia.
   * @param centerRow the row of the center of rotation.
   * @param centerCol the col of the center of rotation.
   * @param rotation  the rotation in degrees.
   * @return the minutia rotated around the given center.
   */
  public static int[] applyRotation(int[] minutia, int centerRow, int centerCol, int rotation) {

	  assert (minutia != null);

	  int x = minutia[1] - centerCol;
	  int y = centerRow - minutia[0];
	  double newX = x * Math.cos(Math.toRadians(rotation)) - y * Math.sin(Math.toRadians(rotation));
	  double newY = x * Math.sin(Math.toRadians(rotation)) + y * Math.cos(Math.toRadians(rotation));
	  int newRow = (int)Math.round(centerRow - newY);
	  int newCol = (int)Math.round(newX + centerCol);
	  int newOrientation = (int)((minutia[2] + rotation) % 360);
	  int[] result = {newRow , newCol, newOrientation};
	  return result;

  }

  /**
   * Applies the specified translation to the minutia.
   *
   * @param minutia        the original minutia.
   * @param rowTranslation the translation along the rows.
   * @param colTranslation the translation along the columns.
   * @return the translated minutia.
   */
  public static int[] applyTranslation(int[] minutia, int rowTranslation, int colTranslation) {

	  assert (minutia != null);

	  int newRow = minutia[0] - rowTranslation;
	  int newCol = minutia[1] - colTranslation;
	  int newOrientation = minutia[2];
	  int[] result = {newRow, newCol, newOrientation};
	  return result;
  } 
  
  /**
   * Computes the row, column, and angle after applying a transformation
   * (translation and rotation).
   *
   * @param minutia        the original minutia.
   * @param centerCol      the column around which the point is rotated.
   * @param centerRow      the row around which the point is rotated.
   * @param rowTranslation the vertical translation.
   * @param colTranslation the horizontal translation.
   * @param rotation       the rotation.
   * @return the transformed minutia.
   */
  public static int[] applyTransformation(int[] minutia, int centerRow, int centerCol, int rowTranslation,
      int colTranslation, int rotation) {

	  assert (minutia != null);

	  int[] newMinutia = applyRotation(minutia, centerRow, centerCol, rotation);
	  newMinutia = applyTranslation(newMinutia, rowTranslation, colTranslation);

	  return newMinutia;
  }

  /**
   * Computes the row, column, and angle after applying a transformation
   * (translation and rotation) for each minutia in the given list.
   *
   * @param minutiae       the list of minutiae.
   * @param centerCol      the column around which the point is rotated.
   * @param centerRow      the row around which the point is rotated.
   * @param rowTranslation the vertical translation.
   * @param colTranslation the horizontal translation.
   * @param rotation       the rotation.
   * @return the list of transformed minutiae.
   */
  public static List<int[]> applyTransformation(List<int[]> minutiae, int centerRow, int centerCol, int rowTranslation,
      int colTranslation, int rotation) {

	  assert (minutiae != null);

	  List<int[]> newMinutiae = new ArrayList<int[]>();

	  for(int m = 0; m < minutiae.size(); ++m){
		  newMinutiae.add(applyTransformation(minutiae.get(m), centerRow, centerCol, rowTranslation, colTranslation, rotation));
	  }
	return newMinutiae;
  }
  /**
   * Counts the number of overlapping minutiae.
   *
   * @param minutiae1      the first set of minutiae.
   * @param minutiae2      the second set of minutiae.
   * @param maxDistance    the maximum distance between two minutiae to consider
   *                       them as overlapping.
   * @param maxOrientation the maximum difference of orientation between two
   *                       minutiae to consider them as overlapping.
   * @return the number of overlapping minutiae.
   */
  public static int matchingMinutiaeCount(List<int[]> minutiae1, List<int[]> minutiae2, int maxDistance,
      int maxOrientation) {

	  assert (minutiae1 != null);
	  assert (minutiae2 != null);

	  int count = 0;

	  for(int indexMinutia1 = 0; indexMinutia1 < minutiae1.size(); ++indexMinutia1 ){
		  int[] minutia1 = minutiae1.get(indexMinutia1);

		  boolean found = false;
		  int indexMinutia2 = 0;

		  while(!found && indexMinutia2 <  minutiae2.size()){

			  int[] minutia2 = minutiae2.get(indexMinutia2);

			  double euclideanDistance = Math.sqrt(Math.pow((minutia1[0] - minutia2[0]), 2) + Math.pow((minutia1[1] - minutia2[1]), 2));

			  if(euclideanDistance <= maxDistance && Math.abs(minutia1[2] - minutia2[2]) <= maxOrientation){
				  ++count;
				  found = true;
			  }
			  ++indexMinutia2;
		  }

	  }

	  return count;
  }

  /**
   * Compares the minutiae from two fingerprints.
   *
   * @param minutiae1 the list of minutiae of the first fingerprint.
   * @param minutiae2 the list of minutiae of the second fingerprint.
   * @return Returns <code>true</code> if they match and <code>false</code>
   *         otherwise.
   */
  public static boolean match(List<int[]> minutiae1, List<int[]> minutiae2) {

	  assert (minutiae1 != null);
	  assert (minutiae2 != null);

		for(int m1 = 0; m1 < minutiae1.size(); ++m1){
			int centerRow = minutiae1.get(m1)[0];
			int centerCol = minutiae1.get(m1)[1];

			for(int m2 = 0; m2 < minutiae2.size(); ++m2){

				int rowTranslation = minutiae2.get(m2)[0] - minutiae1.get(m1)[0];
				int colTranslation = minutiae2.get(m2)[1] - minutiae1.get(m1)[1];

				int rotation =  Math.abs(minutiae2.get(m2)[2] - minutiae1.get(m1)[2]);

				for(int rota = rotation - MATCH_ANGLE_OFFSET; rota <= rotation + MATCH_ANGLE_OFFSET; ++rota){

					List<int[]> newMinutiae = new ArrayList<int[]>();
					newMinutiae = applyTransformation(minutiae2, centerRow, centerCol, rowTranslation, colTranslation, rota);

					if(matchingMinutiaeCount(minutiae1, newMinutiae, DISTANCE_THRESHOLD, ORIENTATION_THRESHOLD) >= FOUND_THRESHOLD){
						return true;
					}
				}
			}
		}
	  return false;
  }

}
