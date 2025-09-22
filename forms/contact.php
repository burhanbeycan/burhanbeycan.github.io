<?php
/**
 * Contact Form Handler for Burhan Beycan CV Website
 * This is a basic PHP contact form handler
 * For production use, consider using more robust solutions
 */

// Only process POST requests
if ($_SERVER["REQUEST_METHOD"] != "POST") {
    http_response_code(405);
    echo "Method Not Allowed";
    exit;
}

// Get form data
$name = isset($_POST['name']) ? strip_tags(trim($_POST['name'])) : '';
$email = isset($_POST['email']) ? trim($_POST['email']) : '';
$subject = isset($_POST['subject']) ? strip_tags(trim($_POST['subject'])) : '';
$message = isset($_POST['message']) ? strip_tags(trim($_POST['message'])) : '';

// Validate form data
$errors = array();

if (empty($name)) {
    $errors[] = 'Name is required.';
}

if (empty($email)) {
    $errors[] = 'Email is required.';
} elseif (!filter_var($email, FILTER_VALIDATE_EMAIL)) {
    $errors[] = 'Email is not valid.';
}

if (empty($subject)) {
    $errors[] = 'Subject is required.';
}

if (empty($message)) {
    $errors[] = 'Message is required.';
}

// If there are errors, return them
if (!empty($errors)) {
    http_response_code(400);
    echo implode('; ', $errors);
    exit;
}

// Email configuration
$to = 'burhanbeycan@hotmail.com'; // Replace with actual email
$email_subject = 'New Contact Form Submission: ' . $subject;

// Create email content
$email_content = "Name: $name\n";
$email_content .= "Email: $email\n\n";
$email_content .= "Subject: $subject\n\n";
$email_content .= "Message:\n$message\n";

// Email headers
$email_headers = "From: $name <$email>\r\n";
$email_headers .= "Reply-To: $email\r\n";
$email_headers .= "Content-Type: text/plain; charset=UTF-8\r\n";

// Send email
if (mail($to, $email_subject, $email_content, $email_headers)) {
    echo 'OK';
} else {
    http_response_code(500);
    echo 'Failed to send message. Please try again later.';
}

// Optional: Log the submission (for debugging)
$log_entry = date('Y-m-d H:i:s') . " - Contact form submission from: $name ($email)\n";
file_put_contents('contact_log.txt', $log_entry, FILE_APPEND | LOCK_EX);
?>
