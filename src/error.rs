#[derive(Debug)]
pub enum Severity {
    Warning,
    Fatal,
}

#[derive(Debug)]
pub struct Error {
    pub severity: Severity,
    pub message: String,
}

impl Error {
    pub fn warning(msg: impl Into<String>) -> Self {
        Self { severity: Severity::Warning, message: msg.into() }
    }
    pub fn fatal(msg: impl Into<String>) -> Self {
        Self { severity: Severity::Fatal, message: msg.into() }
    }
}
