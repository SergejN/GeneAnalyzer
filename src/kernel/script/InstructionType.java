/*
    File:
        InstructionType.java 
 *   
    Revision:
        1.0.0.1
 * 
    Description:
        Specifies the instruction type.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package kernel.script;


public enum InstructionType 
{
    Load,
    Unload,
    Save,
    Filter,
    Analyze,
    DPGPLoad,
    System
};
