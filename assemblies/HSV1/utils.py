from . import expression


def classify(RNA: str) -> str:
    return expression.stages.get(RNA, "Uncharacterized")


palette = {
    "Immediate Early": "#00A2A8",
    "Early": "#BA3BFC",
    "Late": "#FAA623",
    "Latency": "#BA4040",
    "Uncharacterized": "#B4B4B4"
}
