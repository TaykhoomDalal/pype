from .metadata import add_phenotype_metadata
from .mr import analyze as mendelian_randomization
from .phewas import phenome_wide_association


__version__ = "0.1.0"
__all__ = [
	"add_phenotype_metadata",
	"mendelian_randomization",
	"phenome_wide_association",
]