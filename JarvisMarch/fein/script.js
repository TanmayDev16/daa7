document.addEventListener('DOMContentLoaded', () => {
    const canvas = document.getElementById('myCanvas');
    const ctx = canvas.getContext('2d');
    const points = [];

    canvas.addEventListener('click', (event) => {
        const rect = canvas.getBoundingClientRect();
        const x = event.clientX - rect.left;
        const y = event.clientY - rect.top;
        points.push({ x, y });

        drawPoint(ctx, x, y);
        displayCoordinates(points);
        logCoordinates(points);
    });

    function drawPoint(ctx, x, y) {
        ctx.beginPath();
        ctx.arc(x, y, 5, 0, 2 * Math.PI);
        ctx.fillStyle = 'blue';
        ctx.fill();
        ctx.closePath();
    }

    function displayCoordinates(points) {
        const coordinatesList = document.getElementById('coordinatesList');
        coordinatesList.innerHTML = '';
        points.forEach((point, index) => {
            const li = document.createElement('li');
            li.textContent = `Point ${index + 1}: (${point.x}, ${point.y})`;
            coordinatesList.appendChild(li);
        });
    }

    function logCoordinates(points) {
        console.log('Coordinates:', points);
    }
});

