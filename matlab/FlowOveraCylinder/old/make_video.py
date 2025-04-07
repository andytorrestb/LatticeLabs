import cv2
import os
from natsort import natsorted

def images_to_video(image_folder, output_video_path, fps=30, codec='mp4v'):
    # Get image file paths and sort them naturally (e.g., img1, img2, ..., img10)
    image_files = [f for f in os.listdir(image_folder)
                   if f.lower().endswith(('.png', '.jpg', '.jpeg'))]
    image_files = natsorted(image_files)

    if not image_files:
        raise ValueError("No images found in the specified folder.")

    # Read the first image to get dimensions
    first_image_path = os.path.join(image_folder, image_files[0])
    frame = cv2.imread(first_image_path)
    height, width, _ = frame.shape

    # Define the codec and create VideoWriter object
    fourcc = cv2.VideoWriter_fourcc(*codec)
    out = cv2.VideoWriter(output_video_path, fourcc, fps, (width, height))

    # Process and write each image
    for image_file in image_files:
        img_path = os.path.join(image_folder, image_file)
        frame = cv2.imread(img_path)
        if frame.shape[:2] != (height, width):
            frame = cv2.resize(frame, (width, height))  # Resize to match first frame
        out.write(frame)

    out.release()
    print(f"Video saved as {output_video_path}")

# Example usage
if __name__ == "__main__":
    images_to_video(
        image_folder="results-MLS-2025-03-31T21-04-55",  # Folder containing input images
        output_video_path="output_video_1.mp4",
        fps=12
    )
