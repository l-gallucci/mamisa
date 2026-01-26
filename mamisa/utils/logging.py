"""
Logging and output formatting utilities
"""

import sys
from typing import Optional


def log_info(message: str, prefix: str = "[INFO]"):
    """Print info message to stderr"""
    print(f"{prefix} {message}", file=sys.stderr)


def log_error(message: str, prefix: str = "[ERROR]"):
    """Print error message to stderr"""
    print(f"{prefix} {message}", file=sys.stderr)


def log_warning(message: str, prefix: str = "[WARN]"):
    """Print warning message to stderr"""
    print(f"{prefix} {message}", file=sys.stderr)


def log_debug(message: str, prefix: str = "[DEBUG]"):
    """Print debug message to stderr"""
    print(f"{prefix} {message}", file=sys.stderr)


def print_header(title: str, width: int = 70, char: str = "="):
    """Print a formatted header"""
    print(char * width, file=sys.stderr)
    padding = (width - len(title) - 2) // 2
    print(f"{' ' * padding} {title}", file=sys.stderr)
    print(char * width, file=sys.stderr)


def print_section(title: str, width: int = 70, char: str = "="):
    """Print a section header"""
    print(f"\n{char * width}", file=sys.stderr)
    print(f"  {title}", file=sys.stderr)
    print(f"{char * width}", file=sys.stderr)


def format_number(n: int, width: int = 10) -> str:
    """Format number with thousand separators"""
    return f"{n:>{width},}"


def format_percentage(value: float, total: float, decimals: int = 2) -> str:
    """Calculate and format percentage"""
    if total == 0:
        return "0.00%"
    pct = 100 * value / total
    return f"{pct:.{decimals}f}%"


def format_bytes(n_bytes: int) -> str:
    """Format bytes in human-readable format"""
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if n_bytes < 1024.0:
            return f"{n_bytes:.2f} {unit}"
        n_bytes /= 1024.0
    return f"{n_bytes:.2f} PB"


class ProgressBar:
    """Simple progress bar for terminal"""
    
    def __init__(self, total: int, prefix: str = "", width: int = 50):
        self.total = total
        self.current = 0
        self.prefix = prefix
        self.width = width
    
    def update(self, n: int = 1):
        """Update progress by n steps"""
        self.current += n
        self._display()
    
    def _display(self):
        """Display current progress"""
        if self.total == 0:
            return
        
        pct = self.current / self.total
        filled = int(self.width * pct)
        bar = '█' * filled + '░' * (self.width - filled)
        
        print(f"\r{self.prefix} |{bar}| {pct*100:.1f}% ({self.current}/{self.total})", 
              end='', file=sys.stderr)
        
        if self.current >= self.total:
            print(file=sys.stderr)
    
    def finish(self):
        """Mark progress as complete"""
        self.current = self.total
        self._display()
