"""
LLM Reasoning Module for MMSage-Insight
Provides mechanism reasoning using OpenAI API.

Follows rep1212 pattern for expert system LLM integration:
- Structured context building
- Prompt engineering for scientific reasoning
- Response parsing and formatting
"""
import os
import json
from typing import Optional, Dict, List, Any
from dataclasses import dataclass, field

try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False


# Configuration — set OPENAI_API_KEY (and optionally OPENAI_API_BASE, OPENAI_MODEL)
OPENAI_API_KEY = os.getenv("OPENAI_API_KEY", "") or os.getenv("AI2API_KEY", "")
OPENAI_MODEL = os.getenv("OPENAI_MODEL", "gpt-4o")
_default_base = "https://api.openai.com/v1" if OPENAI_API_KEY else ""
OPENAI_API_BASE = os.getenv("OPENAI_API_BASE", "").strip() or _default_base


@dataclass
class ReasoningResult:
    """Result from LLM mechanism reasoning."""
    hypothesis: str
    mechanism: str
    key_insights: List[Dict[str, str]]
    confidence_assessment: str
    limitations: List[str]
    generated_by: str  # 'openai', 'template'
    model: Optional[str] = None
    raw_response: Optional[str] = None
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            'hypothesis': self.hypothesis,
            'mechanism': self.mechanism,
            'key_insights': self.key_insights,
            'confidence_assessment': self.confidence_assessment,
            'limitations': self.limitations,
            'generated_by': self.generated_by,
            'model': self.model
        }


class LLMReasoning:
    """
    LLM-based mechanism reasoning using OpenAI API.
    """
    
    SYSTEM_PROMPT = """You are an expert microbiome-metabolomics researcher specializing in gut microbiota and their metabolic interactions. Your role is to analyze microbe-metabolite associations and propose mechanistic hypotheses.

When analyzing, you should:
1. Consider the biological plausibility of the proposed mechanism
2. Reference relevant metabolic pathways
3. Acknowledge uncertainty and limitations
4. Provide clear, structured reasoning

Format your response as JSON with the following structure:
{
    "hypothesis": "One-sentence summary of the proposed causal relationship",
    "mechanism": "2-3 paragraph detailed explanation of the proposed molecular mechanism",
    "key_insights": [
        {"type": "direction", "content": "Description of causal direction"},
        {"type": "pathway", "content": "Key metabolic pathway involved"},
        {"type": "evidence", "content": "Main supporting evidence"}
    ],
    "confidence_assessment": "Assessment of confidence level with justification",
    "limitations": ["Limitation 1", "Limitation 2"]
}"""

    USER_PROMPT_TEMPLATE = """Analyze the mechanistic relationship between the following microbe and metabolite based on the provided evidence.

{context}

Please provide your analysis in the specified JSON format."""

    def __init__(self, api_key: Optional[str] = None, model: Optional[str] = None):
        """
        Initialize LLM reasoning module.
        
        Args:
            api_key: OpenAI API key (uses env var if not provided)
            model: Model name (uses env var or default if not provided)
        """
        self.api_key = api_key or OPENAI_API_KEY
        self.model = model or OPENAI_MODEL
        self._validate_setup()
    
    def _validate_setup(self) -> None:
        """Validate that the module is properly configured."""
        self.is_configured = bool(self.api_key) and HAS_REQUESTS
        if not HAS_REQUESTS:
            print("Warning: requests library not installed. LLM API calls will not work.")
        if not self.api_key:
            print("Warning: OPENAI_API_KEY not set. LLM reasoning is disabled.")
    
    def reason(
        self,
        context: str,
        bacteria: str,
        metabolite: str,
        use_llm: bool = True
    ) -> ReasoningResult:
        """
        Generate mechanism reasoning.
        
        Args:
            context: Formatted knowledge context from KnowledgeFormatter
            bacteria: Target bacteria name
            metabolite: Target metabolite name
            use_llm: Whether to attempt LLM API call
            
        Returns:
            ReasoningResult with hypothesis and mechanism
        """
        if not use_llm:
            raise RuntimeError("LLM reasoning disabled. Set use_llm=true to run inference.")

        if not self.is_configured:
            raise RuntimeError("OpenAI API is not configured. Set OPENAI_API_KEY to enable LLM reasoning.")

        return self._call_openai(context, bacteria, metabolite)
    
    def _call_openai(
        self,
        context: str,
        bacteria: str,
        metabolite: str
    ) -> ReasoningResult:
        """Call OpenAI API for mechanism reasoning."""
        if not HAS_REQUESTS:
            raise RuntimeError("requests library not available")
        
        headers = {
            "Content-Type": "application/json",
            "Authorization": f"Bearer {self.api_key}"
        }
        
        user_message = self.USER_PROMPT_TEMPLATE.format(context=context)
        
        payload = {
            "model": self.model,
            "messages": [
                {"role": "system", "content": self.SYSTEM_PROMPT},
                {"role": "user", "content": user_message}
            ],
            "temperature": 0.7,
            "max_tokens": 2000,
            "response_format": {"type": "json_object"}
        }
        
        response = requests.post(
            f"{OPENAI_API_BASE}/chat/completions",
            headers=headers,
            json=payload,
            timeout=60
        )
        
        if response.status_code != 200:
            raise RuntimeError(f"OpenAI API error: {response.status_code} - {response.text}")
        
        result = response.json()
        content = result['choices'][0]['message']['content']
        
        # Parse JSON response
        parsed = json.loads(content)
        
        return ReasoningResult(
            hypothesis=parsed.get('hypothesis', 'No hypothesis generated'),
            mechanism=parsed.get('mechanism', 'No mechanism generated'),
            key_insights=parsed.get('key_insights', []),
            confidence_assessment=parsed.get('confidence_assessment', 'Unknown'),
            limitations=parsed.get('limitations', []),
            generated_by='openai',
            model=self.model,
            raw_response=content
        )
    
    
    def is_api_available(self) -> bool:
        """Check if OpenAI API is configured and available."""
        return self.is_configured
    
    def get_config(self) -> Dict[str, Any]:
        """Get current configuration (without exposing full API key)."""
        return {
            'api_configured': bool(self.api_key),
            'api_key_preview': f"{self.api_key[:8]}..." if self.api_key else None,
            'model': self.model,
            'api_base': OPENAI_API_BASE,
            'requests_available': HAS_REQUESTS
        }


# Singleton instance
_llm: Optional[LLMReasoning] = None


def get_llm_reasoning(api_key: Optional[str] = None) -> LLMReasoning:
    """Get or create LLMReasoning singleton."""
    global _llm
    if _llm is None:
        _llm = LLMReasoning(api_key)
    return _llm


def configure_llm(api_key: str, model: Optional[str] = None) -> LLMReasoning:
    """Configure LLM with new API key."""
    global _llm
    _llm = LLMReasoning(api_key, model)
    return _llm


if __name__ == "__main__":
    # Print config only
    llm = LLMReasoning()
    print("=== LLM Configuration ===")
    print(json.dumps(llm.get_config(), indent=2))
