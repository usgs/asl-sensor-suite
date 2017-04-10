package asl.sensor.input;

/**
 * Enum for specifiying unit types, as in values parsed from a response file
 * @author akearns
 *
 */
public enum Unit {

  /**
   * Specifies that a given piece of data has units of m/s
   */
  VELOCITY(1),
  /**
   * Specifies that a given piece of data has units of m/(s^2)
   */
  ACCELERATION(2),
  /**
   * Specifies that a given piece of data has units of m (i.e., no time unit)
   */
  DISPLACEMENT(0);
  
  int secondsExponent;
  
  Unit(int secondExponent) {
    this.secondsExponent = secondExponent;
  }
  
  /**
   * Returns the number of times that data must be differentiated
   * in order to convert from an inputted unit to this unit (i.e., if the 
   * passed parameter is VELOCITY and this is ACCELERATION, this returns 1, 
   * and if the passed parameter is DISPLACEMENT this returns 2).
   * If the number is negative, then it's the number of times the data must be
   * integrated instead
   * @param convertFrom Unit type to be converted from
   * @return Number of times data must be integrated/differentiated as (signed)
   * integer
   */
  public int getDifferentiations(Unit convertFrom) {
    return this.secondsExponent - convertFrom.getSecondsUnit();
  }
  
  /**
   * Returns the power of seconds in the denominator (i.e., 2 for acceleration,
   * which has units of seconds-squared in the denominator)
   * @return The exponent for the seconds of this unit
   */
  public int getSecondsUnit() {
    return this.secondsExponent;
  }
  
}
