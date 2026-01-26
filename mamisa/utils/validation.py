"""
Input validation utilities
"""

import sys
from pathlib import Path
from typing import Optional


class ValidationError(Exception):
    """Custom exception for validation errors"""
    pass


def validate_file_exists(filepath: Path, description: str = "File") -> Path:
    """
    Validate that a file exists and is readable
    
    Args:
        filepath: Path to file
        description: Description of the file for error messages
        
    Returns:
        The validated Path object
        
    Raises:
        ValidationError: If file doesn't exist or is not readable
    """
    if filepath is None:
        raise ValidationError(f"{description} not specified")
    
    filepath = Path(filepath)
    
    if not filepath.exists():
        raise ValidationError(f"{description} not found: {filepath}")
    
    if not filepath.is_file():
        raise ValidationError(f"{description} is not a file: {filepath}")
    
    if not filepath.stat().st_size > 0:
        raise ValidationError(f"{description} is empty: {filepath}")
    
    return filepath


def validate_dir_exists(dirpath: Path, description: str = "Directory") -> Path:
    """
    Validate that a directory exists and is accessible
    
    Args:
        dirpath: Path to directory
        description: Description of the directory for error messages
        
    Returns:
        The validated Path object
        
    Raises:
        ValidationError: If directory doesn't exist or is not accessible
    """
    if dirpath is None:
        raise ValidationError(f"{description} not specified")
    
    dirpath = Path(dirpath)
    
    if not dirpath.exists():
        raise ValidationError(f"{description} not found: {dirpath}")
    
    if not dirpath.is_dir():
        raise ValidationError(f"{description} is not a directory: {dirpath}")
    
    return dirpath


def validate_output_path(filepath: Path, overwrite: bool = False) -> Path:
    """
    Validate output file path
    
    Args:
        filepath: Path to output file
        overwrite: Whether to allow overwriting existing files
        
    Returns:
        The validated Path object
        
    Raises:
        ValidationError: If path is invalid or file exists and overwrite is False
    """
    if filepath is None:
        raise ValidationError("Output path not specified")
    
    filepath = Path(filepath)
    
    # Check if parent directory exists
    if not filepath.parent.exists():
        raise ValidationError(f"Output directory does not exist: {filepath.parent}")
    
    # Check if file already exists
    if filepath.exists() and not overwrite:
        raise ValidationError(f"Output file already exists: {filepath}\n"
                            "Use --force to overwrite")
    
    return filepath


def validate_positive_int(value: int, name: str = "Value", 
                         min_value: int = 1) -> int:
    """
    Validate that a value is a positive integer
    
    Args:
        value: Value to validate
        name: Name of the value for error messages
        min_value: Minimum allowed value (default: 1)
        
    Returns:
        The validated value
        
    Raises:
        ValidationError: If value is not positive
    """
    if not isinstance(value, int):
        raise ValidationError(f"{name} must be an integer")
    
    if value < min_value:
        raise ValidationError(f"{name} must be >= {min_value}, got {value}")
    
    return value


def validate_percentage(value: float, name: str = "Value") -> float:
    """
    Validate that a value is a valid percentage (0-100)
    
    Args:
        value: Value to validate
        name: Name of the value for error messages
        
    Returns:
        The validated value
        
    Raises:
        ValidationError: If value is not in range 0-100
    """
    if not isinstance(value, (int, float)):
        raise ValidationError(f"{name} must be a number")
    
    if not 0 <= value <= 100:
        raise ValidationError(f"{name} must be between 0 and 100, got {value}")
    
    return float(value)


def validate_choice(value: str, choices: list, name: str = "Value") -> str:
    """
    Validate that a value is in a list of allowed choices
    
    Args:
        value: Value to validate
        choices: List of allowed values
        name: Name of the value for error messages
        
    Returns:
        The validated value
        
    Raises:
        ValidationError: If value is not in choices
    """
    if value not in choices:
        raise ValidationError(
            f"{name} must be one of {choices}, got '{value}'"
        )
    
    return value


def check_dependencies(executables: list) -> dict:
    """
    Check if required executables are available in PATH
    
    Args:
        executables: List of executable names to check
        
    Returns:
        Dictionary mapping executable names to their paths (or None if not found)
    """
    import shutil
    
    results = {}
    for exe in executables:
        path = shutil.which(exe)
        results[exe] = path
    
    return results


def validate_fasta_format(filepath: Path) -> bool:
    """
    Basic validation that a file appears to be in FASTA format
    
    Args:
        filepath: Path to file to validate
        
    Returns:
        True if file appears to be FASTA format
        
    Raises:
        ValidationError: If file is not valid FASTA
    """
    from .fasta import open_file
    
    has_header = False
    has_sequence = False
    
    with open_file(filepath, 'rt') as f:
        for i, line in enumerate(f):
            line = line.strip()
            
            if i == 0 and not line.startswith('>'):
                raise ValidationError(
                    f"File does not appear to be FASTA format: {filepath}\n"
                    "First line should start with '>'"
                )
            
            if line.startswith('>'):
                has_header = True
            elif line and has_header:
                has_sequence = True
                break
            
            if i > 100:  # Check first 100 lines only
                break
    
    if not has_header or not has_sequence:
        raise ValidationError(
            f"File does not appear to be valid FASTA format: {filepath}"
        )
    
    return True
